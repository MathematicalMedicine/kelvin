test :
	- rm -rf *.out
	@echo -e "\n\
    **** THE FOLLOWING TEXT SHOULD APPEAR AND IS NOT AN ERROR ****\n\
    make: [test] Error 1 (ignored)\n\
    grep "Some pedigrees have disconnected individuals" kelvin.out\n\
    ERROR - EXITING (read_pedfile.c:225), Some pedigrees have disconnected individuals\n"
	- $(TEST_KELVIN) kelvin.conf >kelvin.out 2>&1
	grep "Some pedigrees have disconnected individuals" kelvin.out

test-USE_DL :

clean :
	- rm -rf *.out *.so *memory.dat compiled
