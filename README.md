# dxMAF

**DNAnexus platform apps for working with MAF cross-species genome alignments**

On the platform, imported MAF alignments reside the public project [MAF whole-genome alignments](https://platform.dnanexus.com/projects/B4Kb8P99KQQ1k0F3Qj2002k9/data/), and [maf_stitcher](https://platform.dnanexus.com/app/maf_stitcher) is a published app. There is also a public project with [stitched alignments of human CCDS ORFs](https://platform.dnanexus.com/projects/B67F658fZ7JX4k0vJxjQ01XZ/data), along with the workflow to generate them.

Build dependencies: [dx-toolkit](http://wiki.dnanexus.com/Downloads), [dx-ocaml](https://github.com/mlin/dx-ocaml)

Applets:

* ``maf_importer``: import .maf[.gz] files into a GTable with a genomic range index
* ``maf_stitcher``: stitch contiguous alignments of given genomic annotations/intervals, producing a GTable with type ``CrossSpeciesAlignments``
* ``CrossSpeciesAlignments_Longest_Isoforms``: given a ``CrossSpeciesAlignments`` table that may contain several overlapping isoforms of each "gene", produce a similar table containing only the lengthiest isoform

These applets use the OCaml library found in ``lib/`` which must be built and installed first.

