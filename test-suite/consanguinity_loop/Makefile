test :
	- rm -f br.out ppl.out kelvin.out
	- time ../../kelvin kelvin.conf > kelvin.out 2>&1
	grep "Not all loops in pedigrees are broken" kelvin.out

clean :
	- rm -f *.out *memory.dat *~
