../resources/importMAF: _tags importMAF.ml
	echo "OCAMLPATH=$(OCAMLPATH)"
	ocamlbuild -use-ocamlfind importMAF.native
	mkdir -p ../resources
	cp -f importMAF.native ../resources/importMAF

clean:
	ocamlbuild -clean
	rm -f ../resources/importMAF
