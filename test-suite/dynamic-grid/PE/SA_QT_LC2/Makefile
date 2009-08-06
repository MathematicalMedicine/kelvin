test :
	- rm -rf *.out
	(time $(TEST_KELVIN) kelvin.conf) >kelvin.out 2>&1
	perl $(KELVIN_ROOT)/utils/wordDiff.pl br.out br.out-baseline .000001

test-FIXED :
	echo Fixed-grid comparison not available for LC runs yet

clean :
	- rm -rf *.out *.so *memory.dat compiled LRTest*

