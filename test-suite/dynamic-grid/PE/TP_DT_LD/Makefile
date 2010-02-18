test :
	- rm -rf *.out
	(time $(TEST_KELVIN) kelvin.conf) >kelvin.out 2>&1
	diff -i br.out br.out-baseline
	diff -i ppl.out ppl.out-baseline

test-USE_DL :
	- rm -rf *.out *.so compiled
	(time $(TEST_KELVIN)-POLYCOMP_DL kelvin.conf) >kelvin-COMP_DL.out 2>&1
	diff -i br.out br.out-baseline
	diff -i ppl.out ppl.out-baseline
	(time $(TEST_KELVIN)-POLYUSE_DL kelvin.conf) >kelvin-USE_DL.out 2>&1
	diff -i br.out br.out-baseline
	diff -i ppl.out ppl.out-baseline

test-FIXED :
	-rm -f LRTest*
	(perl $(KELVIN_ROOT)/LRTest.pl) >kelvin-FIXED.out 2>&1

clean :
	- rm -rf *.out *.so *memory.dat compiled LRTest*
