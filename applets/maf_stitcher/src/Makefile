../resources/maf_stitcher: _tags maf_stitcher.ml
	ocamlbuild -use-ocamlfind maf_stitcher.native
	mkdir -p ../resources
	cp -f maf_stitcher.native ../resources/maf_stitcher

clean:
	ocamlbuild -clean
	rm -f ../resources/maf_stitcher
