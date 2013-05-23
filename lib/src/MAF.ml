open Batteries

include MAFDefs

(* break the input down into "paragraphs" separated by blank lines, stripping out any comments in between. *)
let paragraphs input =
	let lines = IO.lines_of input
	Enum.from
		fun () ->
			let rec skip_blank () =
				match Enum.peek lines with
					| None -> raise Enum.No_more_elements
					| Some line ->
						if (String.length (String.trim line) = 0 ||
							(String.length line > 0 && line.[0] = '#')) then
							Enum.junk lines
							skip_blank ()
			skip_blank ()
			let para = Buffer.create 256
			let rec input_block () =
				match Enum.get lines with
					| None -> ()
					| Some line ->
						if String.length line > 0 then
							Buffer.add_string para line
							Buffer.add_char para '\n'
							input_block ()
			input_block ()
			if Buffer.length para = 0 then raise Enum.No_more_elements
			Buffer.contents para

let block_of_paragraph ?(minimal=false) txt =
	let block = (if minimal then MAFBlockParser.parse_minimally
					else MAFBlockParser.parse) MAFBlockLexer.token (Lexing.from_string txt)
	assert (Array.length block.sequences > 0)
	let len = String.length block.sequences.(0).text
	if block.sequences |> Array.enum |> exists (fun { text } -> String.length text <> len) then
		failwith "inconsistent sequence lengths"
	{ block with unparsed = txt }

let parse input = Enum.map block_of_paragraph (paragraphs input)

let parse_minimally input = Enum.map (block_of_paragraph ~minimal:true) (paragraphs input)

exception Discontiguity of string*int
let stitch ?(src_transform=identity) ?missing_char ?order ~src ~lo ~hi blocks =
	if hi < lo then invalid_arg "Alignment.MAF.stitch: hi < lo"

	let to_splice = ref []
	let pos = ref lo

	List.iter
		fun block ->
			let block_lo = block.sequences.(0).start
			let breflen = block.sequences.(0).size
			let block_hi = block_lo + breflen - 1

			(* Printf.printf "%d %d %d %d %d %s\n" lo block_lo block_hi reflen (Array.length block.sequences) block.sequences.(0).src *)

			if block_hi >= !pos && block_lo <= hi then
				let ids = Array.map src_transform (sources block)
				if ids.(0) <> src_transform src then
					invalid_arg (Printf.sprintf "Alignment.MAF.stitch: block with wrong reference src (expected '%s', got '%s')" ids.(0) (src_transform src))
				if block_lo > !pos then
					raise (Discontiguity (ids.(0), block_lo))
				else if !pos > lo && block_lo < !pos then
					failwith "Alignment.MAF.stitch: detected overlap"
				let conseqs = MAFDefs.to_aln block
				let brefseq,(`RefToCol ref2con),_ = Aln.refseq conseqs
				if String.length brefseq <> breflen then
					failwith "Alignment.MAF.stitch: mismatch between block sequence size and actual reference length"

				(* truncate conseqs to beginning/end of desired region, if necessary *)
				let desired_lo = !pos - block_lo
				let desired_hi = min (hi-block_lo) (breflen-1)
				let desired_length = desired_hi-desired_lo+1
				(*Printf.printf "%d %d %d\n" desired_lo desired_hi desired_length*)
				let conseqs' =
					if desired_lo > lo || desired_hi < block_hi then
						let desired_conlo = if desired_lo > 0 then ref2con.(desired_lo) else 0
						let desired_conhi =
							if desired_hi < breflen-1 then
								ref2con.(desired_hi)
							else
								String.length conseqs.(0) - 1
						Aln.sub_cols conseqs desired_conlo (desired_conhi - desired_conlo + 1)
					else
						conseqs

				to_splice := (ids,conseqs') :: !to_splice
				pos := block_lo+desired_lo+desired_length
		List.sort (fun b1 b2 -> compare b1.sequences.(0).start b2.sequences.(0).start) blocks

	if !pos-1 <> hi then
		raise (Discontiguity ("", !pos))

	Aln.stitch ?missing_char ?order (List.rev !to_splice)
