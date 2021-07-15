all:
	 g++ -fpic -O3 -I/usr/local/igraph -L/usr/local/lib -ligraph main.cc bdgecc.cc -o bound_ecc

debug:
	 g++ -fpic -g3 -I/usr/local/igraph -L/usr/local/lib -ligraph main.cc bdgecc.cc -o bound_ecc
clean:
	rm ./bound_ecc
