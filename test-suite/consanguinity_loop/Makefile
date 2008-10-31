test :
	- rm -f br.out ppl.out kelvin.out
	@echo -e "\n\
    **** THE FOLLOWING TEXT SHOULD APPEAR AND IS NOT AN ERROR ****\n\
    make: [test] Error 255 (ignored)\n\
    grep "Not all loops in pedigrees are broken" kelvin.out\n\
    kelvin.c (653): Not all loops in pedigrees are broken.\n"
	- time ../../kelvin kelvin.conf >kelvin.out 2>&1
	grep "Not all loops in pedigrees are broken" kelvin.out

clean :
	- rm -f *.so *.out *memory.dat *~
