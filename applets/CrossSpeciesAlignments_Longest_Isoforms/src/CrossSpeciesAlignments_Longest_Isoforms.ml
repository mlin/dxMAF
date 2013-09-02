open Batteries
open Printf
open JSON.Operators
open DNAnexus
open Genomics

Printexc.record_backtrace true

(* abbreviated JSON casts *)
module J = struct
  include JSON
  let of_str_list lst = of_list (List.map (fun x -> `String x) lst)
  let of_str_assoc lst = of_assoc (List.map (fun (k,v) -> k, `String v) lst)

exception TransSplicing of string

module Alignment = struct
  type t = (string*int*int*int*bool*GTable.row)
  let lo (_,x,_,_,_,_) = x
  let hi (_,_,x,_,_,_) = x
  let exon_len (_,_,_,x,_,_) = x

  let re_locus = Str.regexp "\\([^:]+\\):\\([0-9]+\\)-\\([0-9]+\\)(\\(.\\))"
  let of_row row =
    (* parse chrom, lo, hi, row out of row. throw on trans-splicing *)
    let id = (match row.(1) with `String s -> s | _ -> assert false)
    let loci =
      (match row.(2) with `String s -> s | _ -> assert false) |> String.nsplit ~by:";" |> List.map
        fun s ->
          Scanf.sscanf s "%s@:%u-%u(%c)"
            fun chrom lo hi strand ->
              assert (lo <= hi)
              chrom, lo, hi, (match strand with '+' -> true | '-' -> false | _ -> assert false)
    let chrom = match List.hd loci with (c,_,_,_) -> c
    loci |> List.iter
      function
        | (c,_,_,_) when c <> chrom -> raise (TransSplicing id)
        | _ -> ()
    let lo = loci |> List.map (function (_,l,_,_) -> l) |> List.fold_left min max_int
    let hi = loci |> List.map (function (_,_,h,_) -> h) |> List.fold_left max min_int
    assert (hi >= lo)
    let strand = match List.hd loci with (_,_,_,s) -> s
    loci |> List.iter
      function
        | (_,_,_,s) when s <> strand -> raise (TransSplicing id)
        | _ -> ()
    let exon_len = loci |> List.map (function (_,l,h,_) -> h-l) |> List.fold_left (+) 0
    (chrom, lo, (hi-1), exon_len, strand, row)

module AlignmentI = IntervalOps.Make(Alignment)

(* main job: instantiate the job tree *)
let main input =
  let alignments_table = GTable.bind_link (input$"alignments")
  GTable.reconfigure ~pagination:1000 alignments_table
  let alignments_desc = GTable.describe alignments_table

  (* download all alignments, load into data structure (print warning for any trans-splicing detected) *)
  let alignments =
    GTable.iterate_rows alignments_table /@ Alignment.of_row |> Enum.fold
      fun m aln ->
        let chrom = match aln with (c,_,_,_,_,_) -> c
        let strand = match aln with (_,_,_,_,s,_) -> s
        let k = (chrom,strand)
        Map.add k (aln :: (try Map.find k m with Not_found -> [])) m
      Map.empty

  (* cluster them*)
  let alignments = Map.map (AlignmentI.cluster % (List.sort AlignmentI.compare)) alignments

  (* filter them *)
  let alignments = 
    alignments |> Map.mapi
      fun (chrom,strand) clusters ->
        let tot = clusters |> List.map List.length |> List.fold_left (+) 0
        printf "%s %c %d %d\n" chrom (if strand then '+' else '-') tot (List.length clusters)
        flush stdout
        clusters |> List.map
          fun cluster ->
            List.fold_left
              fun aln1 aln2 ->
                let l1 = Alignment.exon_len aln1
                let l2 = Alignment.exon_len aln2
                if l1 > l2 then
                  aln1
                else if l1 < l2 then
                  aln2
                else if Random.bool () then
                  aln1
                else 
                  aln2
              List.hd cluster
              List.tl cluster

  (* initialize output GTable *)
  let output_name =
    if input$?"output_name" then J.string (input$"output_name")
    else J.string (alignments_desc$"name") ^ " longest isoforms"
  let output_gtable =
    GTable.make_new
      with_workspace_id
        JSON.of_assoc [
          "name", `String output_name;
          "initializeFrom", J.of_str_assoc [
            "project", JSON.string (alignments_desc$"project");
            "id", GTable.id alignments_table
          ]
        ]
  GTable.reconfigure ~pagination:1000 output_gtable

  (* load filtered alignments into output GTable *)
  foreach (Map.enum alignments)
    fun (_, lst) ->
      foreach (List.enum lst)
        fun (_,_,_,_,_,row) -> GTable.add_row output_gtable (Array.sub row 1 (Array.length row - 1))
  GTable.flush_rows output_gtable

  (* close output GTable *)
  GTable.close output_gtable

  J.of_assoc [
    "longest_alignments", make_link (GTable.id output_gtable)
  ]

(* entry point *)
job_main main
