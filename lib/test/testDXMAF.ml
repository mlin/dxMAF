open Batteries
open Printf
open JSON.Operators

let gt = DNAnexus.GTable.bind (None,"gtable-B4PB620X8Pv0GqQ1Ybp00F2b")

let blocks = DXMAF.gri_query gt ~src:"chr17" ~lo:41_242_946 ~hi:41_243_995

foreach blocks (fun { MAF.unparsed } -> printf "%s\n" unparsed; flush stdout)

let sps_json = JSON.from_string "[
    [\"Human\",\"hg19\"],
    [\"Chimp\",\"panTro2\"],
    [\"Gorilla\",\"gorGor1\"],
    [\"Orangutan\",\"ponAbe2\"],
    [\"Rhesus\",\"rheMac2\"],
    [\"Baboon\",\"papHam1\"],
    [\"Marmoset\",\"calJac1\"],
    [\"Tarsier\",\"tarSyr1\"],
    [\"Mouse_lemur\",\"micMur1\"],
    [\"Bushbaby\",\"otoGar1\"],
    [\"TreeShrew\",\"tupBel1\"],
    [\"Mouse\",\"mm9\"],
    [\"Rat\",\"rn4\"],
    [\"Kangaroo_rat\",\"dipOrd1\"],
    [\"Guinea_Pig\",\"cavPor3\"],
    [\"Squirrel\",\"speTri1\"],
    [\"Rabbit\",\"oryCun2\"],
    [\"Pika\",\"ochPri2\"],
    [\"Alpaca\",\"vicPac1\"],
    [\"Dolphin\",\"turTru1\"],
    [\"Cow\",\"bosTau4\"],
    [\"Horse\",\"equCab2\"],
    [\"Cat\",\"felCat3\"],
    [\"Dog\",\"canFam2\"],
    [\"Microbat\",\"myoLuc1\"],
    [\"Megabat\",\"pteVam1\"],
    [\"Hedgehog\",\"eriEur1\"],
    [\"Shrew\",\"sorAra1\"],
    [\"Elephant\",\"loxAfr3\"],
    [\"Rock_hyrax\",\"proCap1\"],
    [\"Tenrec\",\"echTel1\"],
    [\"Armadillo\",\"dasNov2\"],
    [\"Sloth\",\"choHof1\"],
    [\"Opossum\",\"monDom5\"],
    [\"Wallaby\",\"macEug1\"],
    [\"Platypus\",\"ornAna1\"],
    [\"Chicken\",\"galGal3\"],
    [\"Zebra_finch\",\"taeGut1\"],
    [\"Lizard\",\"anoCar1\"],
    [\"X_tropicalis\",\"xenTro2\"]
  ]"

let sps = sps_json |> JSON.array |> Vect.enum |> List.of_enum |> List.map ((flip ($@)) 0) |> List.map JSON.string |> Array.of_list

let sps_aliases =
  List.fold_left
    fun m entry ->
      match entry |> JSON.array |> Vect.enum |> List.of_enum |> List.map JSON.string with
        | [sp] -> m
        | sp :: aliases -> List.fold_left (fun m alias -> Map.add alias sp m) m aliases
        | _ -> assert false
    Map.empty
    sps_json |> JSON.array |> Vect.enum |> List.of_enum

let src_transform s =
  let asmbl = fst (String.split s ".")
  try Map.find asmbl sps_aliases
  with Not_found -> asmbl

let sps, conseqs = DXMAF.gri_stitch gt ~src_transform ~order:sps ~src:"chr17" ~src_prefix:"hg19." ~lo:41_242_946 ~hi:41_243_995

conseqs |> Array.iteri (fun i seq -> printf ">%s\n%s\n" sps.(i) seq)
