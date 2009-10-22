test :
	- rm -rf *.out
	(time $(TEST_KELVIN) kelvin.conf) >kelvin.out 2>&1
	perl $(KELVIN_ROOT)/utils/wordDiff.pl br.out br.out-baseline .000001
	diff ppl.out ppl.out-baseline

test-USE_DL :
	- rm -rf *.out *.so compiled
	(time $(KELVIN_ROOT)/kelvin-POLYCOMP_DL kelvin.conf) >kelvin-COMP_DL.out 2>&1
	perl $(KELVIN_ROOT)/utils/wordDiff.pl br.out br.out-baseline .000001
	diff ppl.out ppl.out-baseline
	(time $(KELVIN_ROOT)/kelvin-POLYUSE_DL kelvin.conf) >kelvin-USE_DL.out 2>&1
	perl $(KELVIN_ROOT)/utils/wordDiff.pl br.out br.out-baseline .000001
	diff ppl.out ppl.out-baseline

clean :
	- rm -rf *.out *.so *memory.dat compiled
