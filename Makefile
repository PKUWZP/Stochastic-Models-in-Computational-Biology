gillesp_nfkb_sneppen_gene.o: 
	g++ nfkb_decoy_strip.cpp -O3 -o nfkb_decoy_strip.o

clean:
	rm -rf *o nfkb_decoy_strip 
