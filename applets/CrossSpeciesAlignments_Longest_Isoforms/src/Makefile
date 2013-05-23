all:
	mkdir -p ../resources
	$(MAKE) ../resources/CrossSpeciesAlignments_Longest_Isoforms

../resources/CrossSpeciesAlignments_Longest_Isoforms: _tags CrossSpeciesAlignments_Longest_Isoforms.ml
	ocamlbuild -use-ocamlfind CrossSpeciesAlignments_Longest_Isoforms.native
	cp -f CrossSpeciesAlignments_Longest_Isoforms.native ../resources/CrossSpeciesAlignments_Longest_Isoforms

clean:
	ocamlbuild -clean
	rm -f ../resources/CrossSpeciesAlignments_Longest_Isoforms
