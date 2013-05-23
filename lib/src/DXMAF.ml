open Batteries
open DNAnexus
open JSON.Operators

let memoize f =
  let t = Hashtbl.create 16
  fun x ->
    try Hashtbl.find t x
    with
      | Not_found ->
          let y = f x
          Hashtbl.add t x y
          y

let cached_describe =
  memoize
    fun dxid -> GTable.describe (GTable.bind (None,dxid))

let gri_query ?limit gt ~src ~lo ~hi =
  let desc = cached_describe (GTable.id gt)
  let types = JSON.array (desc$"types")
  if not (Vect.mem (`String "maf") types) then
    invalid_arg "DXMAF.gri_query: GTable does not have type maf"
  (* perform GRI query *)
  let rows = GTable.query_rows ?limit gt "gri" (`GRI (src,lo,hi,`Overlap))
  (* map results to parse the MAF blocks *)
  rows |> Enum.map
    function
      | [| `Int rowid; `String src; `Int lo; `Int hi; `String raw_block |] ->
          Option.get (Enum.get (MAF.parse (IO.input_string raw_block)))
      | row ->
          let row_json = JSON.to_string (GTable.json_of_row row)
          failwith ("DXMAF.gri_query: couldn't understand GTable row " ^ row_json)

let gri_stitch ?src_transform ?missing_char ?order ?src_prefix gt ~src ~lo ~hi =
  let blocks = List.of_enum (gri_query gt ~src ~lo ~hi)
  let maf_src = match src_prefix with Some pfx -> pfx ^ src | None -> src
  MAF.stitch ?src_transform ?missing_char ?order ~src:maf_src ~lo ~hi:(hi-1) blocks
