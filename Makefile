default: clean keygen sign verify 
keygen:
	gcc -march=native -mtune=native -funroll-loops -O3 -std=c99 -o keygen keygen.c ntt.c packing.c poly.c shake.c
sign:
	gcc -march=native -mtune=native -funroll-loops -O3 -std=c99 -o sign sign.c ntt.c  packing.c poly.c shake.c
verify:
	gcc -march=native -mtune=native -funroll-loops -O3 -o verify verify.c ntt.c packing.c poly.c shake.c
	
clean:
	-rm $(objects) keygen sign verify

