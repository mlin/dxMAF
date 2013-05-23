open Printf
open Batteries
open JSON.Operators
open DNAnexus

Printexc.record_backtrace true

(* abbreviated JSON casts *)
module J = struct
  include JSON
  let of_str_list lst = of_list (List.map (fun x -> `String x) lst)
  let of_str_assoc lst = of_assoc (List.map (fun (k,v) -> k, `String v) lst)

type region = {
  rgn_id : string;
  rgn_loci : (string*int*int*[`Plus|`Minus]) list
}

(* Index a Spans/gri table, in which each desired region is a single
   genomic interval *)
let index_spans tbl revcomp =
  (* figure out if we have names and strands to work with *)
  let allcols = GTable.columns tbl
  let have_names = Array.mem ("name",`String) allcols
  let have_strands = revcomp && Array.mem ("strand",`String) allcols

  (* specify which columns to retrieve*)
  let columns = match have_names, have_strands with
    | true, true -> [|"name";"chr";"lo";"hi";"strand"|]
    | false, true -> [|"chr";"lo";"hi";"strand"|]
    | true, false -> [|"name";"chr";"lo";"hi"|]
    | false, false -> [|"chr";"lo";"hi"|]

  (* derive a region from an individual row, taking into account what
     information we have to work with (above) *)
  let fresh_region_id =
    let next_region_id = ref 0
    fun () ->
      incr next_region_id
      sprintf "%06d" !next_region_id
  let region_of_row = function
    | [| `String name; `String chr; `Int lo; `Int hi; `String strand |] ->
        assert (have_names && have_strands)
        { rgn_id = name; rgn_loci = [chr,lo,hi,(if strand = "-" then `Minus else `Plus)] }
    | [| `String chr; `Int lo; `Int hi; `String strand |] ->
        assert (not have_names && have_strands)
        { rgn_id = fresh_region_id(); rgn_loci = [chr,lo,hi,(if strand = "-" then `Minus else `Plus)] }
    | [| `String name; `String chr; `Int lo; `Int hi |] ->
        assert (have_names && not have_strands)
        { rgn_id = name; rgn_loci = [chr,lo,hi,`Plus] }
    | [| `String chr; `Int lo; `Int hi |] ->
        assert (not have_names && not have_strands)
        { rgn_id = fresh_region_id(); rgn_loci = [chr,lo,hi,`Plus] }
    | _ -> assert false

  (* Apply that to each row of the table and make a list of the results *)
  GTable.iterate_rows ~columns tbl /@ region_of_row |> List.of_enum

(* Index a Genes table, in which each desired region may be a transcript
   spliced from multiple exons *)
let index_transcripts tbl transcript_type exon_type revcomp =
  let columns = [|"type";"span_id";"parent_id";"name";"chr";"lo";"hi";"strand"|]
  let exons = Hashtbl.create 1024
  let transcripts =
    GTable.iterate_rows ~columns tbl |> Enum.fold
      fun transcripts row ->
        match row with
          | [| `String ty; `Int span_id; `Int parent_id; `String name;
               `String chr; `Int lo; `Int hi; `String strand |] ->
              if ty = transcript_type then (span_id,name) :: transcripts
              else if ty = exon_type then
                let locus = (chr,lo,hi,(if strand = "-" then `Minus else `Plus))
                Hashtbl.add exons parent_id locus
                transcripts
              else transcripts
          | _ -> assert false
      []
  transcripts |> List.rev |> List.filter_map
    fun (span_id, name) ->
      let loci = List.sort compare (Hashtbl.find_all exons span_id)
      if loci = [] then
        eprintf "%s: no exons found, skipping!\n" name; flush stderr
        None
      else if List.length (List.unique (List.map (fun (chr,_,_,_) -> chr) loci)) > 1 then
        eprintf "%s: exons found on multiple chromosomes/contigs, skipping!\n" name; flush stderr
        None
      else
        (* If the region is - strand, reverse the order of the exons here -- we'll take
           care of reverse-complementing the individual exon alignments once we have them
           in hand, below. *)
        let strand =
          if List.for_all (function (_,_,_,`Plus) -> true | _ -> false) loci then Some `Plus
          else if List.for_all (function (_,_,_,`Minus) -> true | _ -> false) loci then Some `Minus
          else None
        match strand with
          | None ->
              eprintf "%s: detected mixed strands, skipping!\n" name; flush stderr
              None
          | Some `Minus when not revcomp ->
              let plus_loci = List.map (fun (chr,lo,hi,_) -> (chr,lo,hi,`Plus)) loci
              Some { rgn_id = name; rgn_loci = plus_loci }
          | Some `Minus -> Some { rgn_id = name; rgn_loci = List.rev loci }
          | Some `Plus -> Some { rgn_id = name; rgn_loci = loci }

(* dispatch to index_spans or index_transcripts, based on the type of the
   given gtable *)
let index_regions regions_table transcript_type exon_type revcomp =
  if Vect.mem (`String "Genes") (J.array (GTable.describe regions_table $ "types")) then
    index_transcripts regions_table transcript_type exon_type revcomp
  else index_spans regions_table revcomp

(* in-place Knuth shuffle *)
let shuffle ar =
  let n = Array.length ar
  for i = 0 to n-2 do
    let j = i + Random.int (n-i)
    let tmp = ar.(i)
    ar.(i) <- ar.(j)
    ar.(j) <- tmp

let partition_regions regions regions_per_job =
  let regions = Array.of_list regions
  let n = Array.length regions

  shuffle regions
  
  (* adjust regions_per_job to make distribution more uniform *)
  let jobs = 1 + (Array.length regions - 1) / regions_per_job
  let k = int_of_float (ceil (float n /. float jobs))

  (0 --^ jobs) /@ (fun i -> Array.to_list (Array.sub regions (i*k) (min k (n-i*k)))) |> List.of_enum

(* generate some functions needed for MAF massaging *)
let make_assembly_map assemblies_json =
  let src_transform =
    let tbl =
      Vect.enum (J.array assemblies_json) |> fold
        fun m entry -> Map.add (J.string (entry$@1)) (J.string (entry$@0)) m
        Map.empty
    tbl |> Map.iter (fun k v -> eprintf "%s => %s\n" k v)
    fun asmbl_seq ->
      try
        Map.find (fst (String.split asmbl_seq ".")) tbl
      with Not_found ->
        asmbl_seq

  let order = Vect.enum (J.array assemblies_json) /@ (fun entry -> J.string (entry$@0)) |> Array.of_enum

  order |> Array.iter (fun sp -> eprintf "%s " sp)
  eprintf "\n"

  let ref_asmbl = J.string (assemblies_json$@0$@1)

  src_transform, order, ref_asmbl

let comp = function
  | 'A' -> 'T'
  | 'a' -> 't'
  | 'G' -> 'C'
  | 'g' -> 'c'
  | 'C' -> 'G'
  | 'c' -> 'g'
  | 'T' -> 'A'
  | 't' -> 'a'
  | 'N' -> 'N'
  | 'n' -> 'n'
  | x -> failwith (sprintf "cannot complement %c" x)

let stitch_region ~maf_gtable ~src_transform ~order ~ref_asmbl region =
  let src_prefix = ref_asmbl ^ "."
  try
    (* retrieve alignment of each locus, reverse-complementing each one if indicated *)
    let locus_alns = 
      region.rgn_loci |> List.map
        fun (chr,lo,hi,strand) ->
          let (sps,fwd_aln) = DXMAF.gri_stitch maf_gtable ~src_transform ~order ~src_prefix ~src:chr ~lo ~hi
          if strand <> `Minus then (sps,fwd_aln)
          else (sps, Aln.revmap comp fwd_aln)
    (* stitch them together *)
    let sps, aln = Aln.stitch locus_alns
    assert (sps = order)

    let loci_txt =
      String.join ";" 
        region.rgn_loci |> List.map
          fun (chr,lo,hi,strand) ->
            sprintf "%s:%d-%d(%c)" chr lo hi
              match strand with `Minus -> '-' | _ -> '+'

    Some (region.rgn_id, loci_txt, aln)
  with
    | MAF.Discontiguity (seq,pos) ->
        eprintf "%s: alignment discontiguity at %s:%d, skipping!\n" region.rgn_id seq pos
        None

let make_alignments_gtable name assemblies applet_id input =
  let species_cols =
    (Vect.enum (J.array assemblies)
      /@ (fun entry -> (J.string (entry$@0)), `String)
      |> List.of_enum)

  GTable.make_new (J.of_assoc [
    "name", `String name;
    "types", J.of_str_list ["CrossSpeciesAlignments"];
    "columns", GTable.json_of_columns (Array.of_list (("_name",`String) :: ("_loci",`String) :: species_cols));
    "indices", GTable.json_of_indices ["_name", `Lexicographic ["_name",`Asc,None]];
    "details", J.of_assoc [
      "maf_stitcher", make_link applet_id;
      "maf_stitcher_input", input
    ]
  ])

type subjob_input = {
  alignments_table_id : string;
  species_set_info : JSON.t;
  subjob_regions : region list
}

(* main job: instantiate the job tree *)
let mainjob input =
  (* load inputs from JSON *)
  let species_set = J.string (input$"species_set")
  assert (String.length species_set > 0)
  let regions_table = GTable.bind_link (input$"regions")
  let regions_per_job = JSON.int (input$"transcripts_per_job")
  let revcomp = not (J.bool (input$"no_revcomp"))

  (* resolve species_set *)
  let job_desc = api_call [job_id (); "describe"] J.empty
  let applet_id =
    if job_desc$?"app" then J.string ((api_call [J.string (job_desc$"app"); "describe"] JSON.empty)$"applet")
    else J.string (job_desc$"applet")
  let applet_deets = api_call [applet_id; "getDetails"] J.empty
  let species_set_info = 
    try applet_deets$"species_sets"$species_set
    with J.No_key _ -> raise (AppError (sprintf "Unknown species set %s. Please see the maf_stitcher documentation for available species sets." species_set))

  (* get MAF *)
  let maf = GTable.bind_link (species_set_info$"maf")

  (* ensure regions and MAF refer to the same reference genome *)
  if not (J.bool (input$"no_reference_check")) then
    let regions_ref_link = get_link (GTable.get_details regions_table $ "original_contigset")
    let maf_ref_link = get_link (GTable.get_details maf $ "original_contigset")
    match regions_ref_link, maf_ref_link with
      | (_, id1), (_, id2) when id1 = id2 -> ()
      | _ -> raise (AppError "The species_set and regions table do not refer to the same reference genome. Set no_reference_check=true to override this check.")

  (* read and index entire regions table *)
  let all_regions =
    index_regions
      regions_table
      JSON.string (input$"transcript_type")
      JSON.string (input$"exon_type")
      revcomp

  eprintf "Found %d transcripts\n" (List.length all_regions)

  let chunks = partition_regions all_regions regions_per_job

  eprintf "Partitioned into %d chunks:" (List.length chunks)

  (* initialize GTable and launch jobs *)
  let output_name =
    if input$?"output_name" then J.string (input$"output_name")
    else sprintf "%s %s alignments" (J.string (GTable.describe regions_table$"name")) species_set
  let gtable = make_alignments_gtable output_name (species_set_info$"assemblies") applet_id input
  let subjobs =
    chunks |> List.map
      fun subjob_regions ->
        eprintf "%d " (List.length subjob_regions)
        (* Marshal subjob input and upload it as a platform file *)
        let subjob_input_fn =
          Batteries.File.with_temporary_out
            fun outfile fn ->
              Marshal.output outfile { alignments_table_id = GTable.id gtable;
                                       species_set_info; subjob_regions }
              fn
        let subjob_input_file = File.upload_new subjob_input_fn
        (* launch subjob with this input *)
        let options =
          if JSON.bool (input$"spot") then
            J.of_assoc [
              "systemRequirements", J.of_assoc [
                "subjob", J.of_assoc [
                  "instanceType", `String "dx_s_m1.xlarge"
                ]
              ]
            ]
          else JSON.empty
        new_job ~options "subjob" (input $+ ("subjob",`String (File.id subjob_input_file)))
  eprintf "\n"

  let finishjob_input =
    (input
      $+ ("finishjob", J.of_assoc [
            "alignments", `String (GTable.id gtable);
            "stats", J.of_list (List.map (fun subjob -> make_jobref subjob "stats") subjobs)]))
  let finishjob = new_job "finish" finishjob_input

  J.of_assoc ["alignments", make_jobref finishjob "alignments"]

(* subjob: process some of the regions files *)
let subjob input =
  let subjob_input_file = File.bind (None, (J.string (input$"subjob")))
  File.download subjob_input_file "subjob_input"
  let subjob_input =
    Batteries.File.with_file_in "subjob_input"
      fun infile -> ((Marshal.input infile):subjob_input)
  let src_transform, order, ref_asmbl = make_assembly_map (subjob_input.species_set_info$"assemblies")

  let worker_process region =
    try
      let maf_gtable = GTable.bind_link (subjob_input.species_set_info$"maf")
      stitch_region ~maf_gtable ~src_transform ~order ~ref_asmbl region
    with
      | DNAnexus.AppError msg ->
          raise
            ForkWork.ChildExn [
              "AppError"; region.rgn_id; msg
            ]
      | exn ->
          raise
            ForkWork.ChildExn [
              "AppInternalError";
              region.rgn_id;
              Printexc.to_string exn;
              Printexc.get_backtrace ()
            ]

  let alignments_gtable = GTable.bind (None,subjob_input.alignments_table_id)
  let stats =
    try 
      (* TODO: pipeline stitching and addRows for individual regions *)
      ForkWork.map_list ~maxprocs:(4*ForkWork.ncores()) ~fail_fast:true worker_process subjob_input.subjob_regions |> List.map
        function
          | Some (rgn_id,loci_txt,aln) ->
              assert (Array.length aln = Array.length order)
              GTable.add_row
                alignments_gtable
                Array.concat [
                  [| `String rgn_id; `String loci_txt |];
                  Array.map (fun x -> `String x) aln
                ]
              0
          | None -> 1
    with
      | ForkWork.ChildExn ["AppError"; region_id; msg] ->
          eprintf "AppError on %s: %s\n" region_id msg
          raise (DNAnexus.AppError (sprintf "%s: %s" region_id msg))
      | ForkWork.ChildExn ["AppInternalError"; region_id; msg; bt] ->
          eprintf "AppInternalError on %s: %s\n" region_id msg
          if bt <> "" then eprintf "%s\n" bt
          raise (DNAnexus.AppInternalError (sprintf "%s: %s" region_id msg))
      | _ -> assert false
  GTable.flush_rows alignments_gtable
  J.of_assoc ["stats", J.of_assoc ["failures", `Int (List.fold_left (+) 0 stats)]]

(* finish-up job: collect stats & initiate GTable closure *)
let finishjob input =
  let gtable = GTable.bind (None, JSON.string (input$"finishjob"$"alignments"))

  let stats =
    Enum.fold
      fun failures j -> failures + JSON.int (j$"failures")
      0
      Vect.enum (JSON.array (input$"finishjob"$"stats"))
  let stats_json = J.of_assoc ["failures", `Int stats]

  GTable.set_details gtable (GTable.get_details gtable $+ ("stitch_stats", stats_json))

  GTable.close gtable
  J.of_assoc ["alignments", make_link (GTable.id gtable)]

(* entry point *)
job_main
  fun input ->
    if input $? "subjob" then subjob input
    else if input $? "finishjob" then finishjob input
    else mainjob input
