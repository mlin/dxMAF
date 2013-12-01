open Printf
open Batteries
open JSON.Operators

(* abbreviated JSON casts *)
let j_lst = JSON.of_list
let j_strlst lst = JSON.of_list (List.map (fun x -> `String x) lst)
let j_obj = JSON.of_assoc
let j_strobj lst = j_obj (List.map (fun (k,v) -> k, `String v) lst)

(* type for collecting & combining stats (chromosomes, ranges, block counts) *)
module Stats = struct
  type t = (string,int*int*int) Map.t
  let empty = Map.empty
  let add ?(n=1) src lo hi x =
    try
      let n', lo', hi' = Map.find src x
      Map.add src ((n+n'), min lo lo', max hi hi') x
    with Not_found -> Map.add src (n,lo,hi) x
  let union x y = Map.foldi (fun src (n,lo,hi) y' -> add ~n src lo hi y') x y
  let to_json x =
    Map.foldi
      fun src (n,lo,hi) j ->
        let info = j_obj [
          "lo", `Int lo;
          "hi", `Int hi;
          "blocks", `Int n
        ]
        j $+ (src,info)
      x
      JSON.empty
  let of_json x =
    List.fold_left
      fun stats src ->
        let n = JSON.int (x$src$"blocks")
        let lo = JSON.int (x$src$"lo")
        let hi = JSON.int (x$src$"hi")
        add ~n src lo hi stats
      empty
      JSON.obj_keys x

(* import one MAF file to the GTable *)
let import_MAF_dxfile ref_genome_name dxfile gtable =
  let ref_genome_name_len = String.length ref_genome_name
  let desc = DNAnexus.File.describe dxfile
  let is_gz = Filename.check_suffix (JSON.string (desc$"name")) ".gz"
  with_dispose ~dispose:(fun fn -> try Sys.remove fn with _ -> ())
    fun fn ->
      (* download and begin streaming the MAF file *)
      DNAnexus.File.download dxfile fn
      let infile =
        if is_gz then
          Unix.open_process_in ("zcat " ^ fn)
        else
          File.open_in fn
      infile |> with_dispose ~dispose:IO.close_in
        fun _ ->
          (* iterate over MAF blocks *)
          MAF.parse infile |> Enum.fold
            fun stats {MAF.sequences; unparsed} ->
              let { MAF.src; start; size } = sequences.(0)

              let chr =
                if ref_genome_name_len > 0 then
                  (* check that reference genome matches *)
                  if (String.length src <= ref_genome_name_len+1
                      || String.left src ref_genome_name_len <> ref_genome_name
                      || src.[ref_genome_name_len] <> '.') then
                    raise (DNAnexus.AppError (sprintf "Expected first sequence in MAF block to originate from %s.<<chromosome>>, but found %s. The MAF file and reference genome may be mismatched." ref_genome_name src))
                  (* then strip it off to get the chromosome/contig name*)
                  String.right src (String.length src - ref_genome_name_len - 1)
                else
                  src

              (* add row for this MAF block *)
              DNAnexus.GTable.add_row gtable [|
                `String chr;
                `Int start;
                `Int (start+size);
                `String unparsed
              |]
              (* update stats *)
              Stats.add chr start (start+size) stats
            Stats.empty
    sprintf "Process%d_Thread%d.maf%s" (Unix.getpid ()) (Thread.id (Thread.self ())) (if is_gz then ".gz" else "")

let make_alignments_gtable name =
  DNAnexus.GTable.make_new (j_obj [
    "name", `String name;
    "types", j_strlst ["gri"; "Spans"; "maf"];
    "columns", j_lst [
      j_strobj ["name", "chr"; "type", "string"];
      j_strobj ["name", "lo"; "type", "int32"];
      j_strobj ["name", "hi"; "type", "int32"];
      j_strobj ["name", "maf_block"; "type", "string"]
    ];
    "indices", j_lst [
      j_strobj [
        "name", "gri";
        "type", "genomic";
        "chr", "chr";
        "lo", "lo";
        "hi", "hi"
      ]
    ]
  ])

open DNAnexus

(* main job: instantiate the job tree *)
let mainjob input =
  let files = Vect.enum (JSON.array (input$"maf_files")) /@ (get_link %> File.bind) |> List.of_enum
  let maxjobs = JSON.int (input$"maxjobs")
  let ref_genome_name =
    if input$?"reference_genome_name" then
      JSON.string (input$"reference_genome_name")
    else
      JSON.string (Record.describe (Record.bind (get_link (input$"reference_genome"))) $ "name")

  (* describe all files, since we'll be using their sizes *)
  let file_descs =
    List.fold_left
      fun map dxfile -> Map.add (File.id dxfile) (File.describe dxfile) map
      Map.empty
      files

  (* sort files largest to smallest *)
  let cmpfiles x y =
    compare
      JSON.int (Map.find (File.id y) file_descs $ "size")
      JSON.int (Map.find (File.id x) file_descs $ "size")
  let files = List.sort cmpfiles files

  (* divvy up the files by assigning the next largest file to the subjob with
     the least total file size; use a heap to keep track of the latter *)
  let empty_jobs = Heap.of_list (Array.to_list (Array.init maxjobs (fun _ -> 0, [])))
  let jobheap = 
    files |> List.fold_left
      fun jobheap file ->
        let filesz = JSON.int (Map.find (File.id file) file_descs $ "size")
        let jobsz, jobfiles = Heap.find_min jobheap
        (* replace the job in the heap to include the new file*)
        Heap.insert (Heap.del_min jobheap) ((jobsz+filesz),(file :: jobfiles))
      empty_jobs
  (* filter out empty jobs *)
  let jobs = List.filter ((<>) []) (List.map snd (Heap.to_list jobheap))

  (* initialize GTable and launch jobs *)
  let gtable = make_alignments_gtable (JSON.string (input$"output_name"))
  let subjobs =
    jobs |> List.map
      fun files ->
        let subjob_input = j_obj [
          "alignments", `String (GTable.id gtable);
          "files", j_lst (List.map (File.id %> make_link) files);
          "reference_genome_name", `String ref_genome_name
        ]
        new_job "subjob" (input $+ ("subjob",subjob_input))

  let finishjob_input =
    (input
      $+ ("stats",j_lst (List.map (fun subjob -> make_jobref subjob "stats") subjobs))
      $+ ("finishjob", j_obj [ "alignments", `String (GTable.id gtable) ]))
  (* wanted to nest "stats" inside "finishjob", but the JM does not understand JBORs 
     at that nesting level *)
  let finishjob = new_job "finish" finishjob_input

  j_obj ["alignments", make_jobref finishjob "alignments"]

(* subjob: process some of the MAF files, collect statistics *)
let subjob input =
  let gtable_id = JSON.string (input$"subjob"$"alignments")
  let ref_genome_name = JSON.string (input$"subjob"$"reference_genome_name")

  let worker_process file_link =
    try
      let gtable = GTable.bind (None,gtable_id)
      GTable.reconfigure ~pagination:10000 gtable
      let ans = import_MAF_dxfile ref_genome_name (File.bind file_link) gtable
      GTable.flush_rows gtable
      ans
    with
      | DNAnexus.AppError msg ->
          raise
            ForkWork.ChildExn [
              "AppError"; Option.default "" (fst file_link); snd file_link; msg
            ]
      | Curl.CurlException(_,curlcode,curlmsg) ->
          raise
            ForkWork.ChildExn [
              "AppInternalError";
              Option.default "" (fst file_link);
              snd file_link;
              sprintf "Curl.CurlException(%d,\"%s\")" curlcode curlmsg;
              Printexc.get_backtrace ()
            ]
      | exn ->
          raise
            ForkWork.ChildExn [
              "AppInternalError";
              Option.default "" (fst file_link);
              snd file_link;
              Printexc.to_string exn;
              Printexc.get_backtrace ()
            ]

  let errmsg file_link msg =
    let file = File.bind file_link
    try
      let nm = JSON.string (File.describe file $ "name")
      sprintf "Error while importing %s (%s): %s." nm (File.id file) msg
    with
      | _ -> sprintf "While importing %s: %s" (File.id file) msg

  let file_links = Vect.enum (JSON.array (input$"subjob"$"files")) /@ get_link |> List.of_enum

  let stats =
    try ForkWork.map_list ~fail_fast:true worker_process file_links
    with
      | ForkWork.ChildExn ["AppError"; pid; fid; msg] ->
          eprintf "AppError on %s: %s\n" fid msg
          raise (DNAnexus.AppError (errmsg ((if pid = "" then None else Some pid),fid) msg))
      | ForkWork.ChildExn ["AppInternalError"; pid; fid; msg; bt] ->
          eprintf "AppInternalError on %s: %s\n" fid msg
          if bt <> "" then eprintf "%s\n" bt
          raise (DNAnexus.AppInternalError (errmsg ((if pid = "" then None else Some pid),fid) msg))
      | _ -> assert false

  j_obj ["stats", Stats.to_json (List.fold_left Stats.union Stats.empty stats)]

(* finish-up job: combine and record statistics, initiate GTable closure *)
let finishjob input =
  let gtable = GTable.bind (None, JSON.string (input$"finishjob"$"alignments"))

  let stats =
    Enum.fold
      fun stats j -> Stats.union (Stats.of_json j) stats
      Stats.empty
      Vect.enum (JSON.array (input$"stats"))

  let deets = j_obj [
    "original_files", (input$"maf_files");
    "original_contigset", (input$"reference_genome");
    "maf_stats", Stats.to_json stats
  ]
  GTable.set_details gtable deets

  GTable.close gtable
  j_obj ["alignments", make_link (GTable.id gtable)]

(* entry point *)
job_main
  fun input ->
    if input $? "subjob" then subjob input
    else if input $? "finishjob" then finishjob input
    else mainjob input
