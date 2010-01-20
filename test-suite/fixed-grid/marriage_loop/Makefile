test :
	- rm -rf *.out
	@echo -e "\n\
    **** THE FOLLOWING TEXT SHOULD APPEAR AND IS NOT AN ERROR ****\n\
    make: [test] Error 255 (ignored)\n\
    grep "Not all loops in pedigrees are broken" kelvin.out\n\
    kelvin.c (653): Not all loops in pedigrees are broken.\n"
	- $(TEST_KELVIN) kelvin.conf >kelvin.out 2>&1
	grep "Not all loops have been broken in pedigrees" kelvin.out

test-USE_DL :

clean :
	- rm -rf *.out *.so *memory.dat compiled
