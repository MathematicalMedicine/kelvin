test :
	- rm -f br.out ppl.out kelvin.out
	(time $(TEST_KELVIN) kelvin.conf) >kelvin.out 2>&1
	perl $(KELVIN_ROOT)/utils/wordDiff.pl br.out br.out-baseline .000001
	perl $(KELVIN_ROOT)/utils/wordDiff.pl ppl.out ppl.out-baseline .000001

test-USE_DL :
	echo Compilation not appropriate for this test.

clean :
	- rm -rf *.so *.out *memory.dat *~ compiled
