(** Multiple Alignment Format (MAF)
@see <http://genome.ucsc.edu/FAQ/FAQformat#format5> MAF format (UCSC)
*)

open Batteries

type sequence = {
  src : string;
  start : int;  (** ZERO-based, inclusive coordinate *)
  size : int;
  strand : [`plus | `minus];
  src_size : int;
  text : string;
}

type block = {
  attributes : (string*string) list;
  sequences : sequence array;

  unparsed : string; (** original text of this MAF block *)
}

(** retrieve the [src] field of each sequence *)
val sources : block -> string array

val to_aln : block -> Aln.t

val parse : IO.input -> block Enum.t

(**/**)
val parse_minimally : IO.input -> block Enum.t
(**/**)

exception Discontiguity of string*int

(** [stitch blocks src lo hi] stitches a contiguous alignment of the interval [\[lo,hi\]] on a
reference sequence [src] given the pertinent MAF blocks. Here [lo] and [hi] are zero-based,
inclusive coordinates. The first sequence of each block in [blocks] must match [src]. The blocks
must be non-overlapping.

@return an alignment block for exactly the interval [\[lo,hi\]] of the reference sequence, with the
name (src) of each row.
@param order as in {!Block.stitch}
*)
val stitch : ?src_transform:(string->string) -> ?missing_char:char -> ?order:(string array) -> src:string -> lo:int -> hi:int -> block list -> (string array)*Aln.t
