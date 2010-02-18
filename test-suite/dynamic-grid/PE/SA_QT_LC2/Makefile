test :
	- rm -rf *.out
	(time $(TEST_KELVIN) kelvin.conf) >kelvin.out 2>&1
	perl $(KELVIN_ROOT)/utils/wordDiff.pl br.out br.out-baseline .000001

test-FIXED :
	echo Fixed-grid comparison not available for LC runs yet

test-USE_DL :
	- rm -rf *.out *.so compiled
	(time $(TEST_KELVIN)-POLYCOMP_DL kelvin.conf) >kelvin-COMP_DL.out 2>&1
	perl $(KELVIN_ROOT)/utils/wordDiff.pl br.out br.out-baseline .000001
	(time $(TEST_KELVIN)-POLYUSE_DL kelvin.conf) >kelvin-USE_DL.out 2>&1
	perl $(KELVIN_ROOT)/utils/wordDiff.pl br.out br.out-baseline .000001

clean :
	- rm -rf *.out *.so *memory.dat compiled LRTest*

