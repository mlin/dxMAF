open Batteries

val gri_query : ?limit:int -> DNAnexus.GTable.t -> src:string -> lo:int -> hi:int -> MAF.block Enum.t

val gri_stitch : ?src_transform:(string->string) -> ?missing_char:char -> ?order:(string array) -> ?src_prefix:string -> DNAnexus.GTable.t -> src:string -> lo:int -> hi:int -> (string array)*Aln.t

