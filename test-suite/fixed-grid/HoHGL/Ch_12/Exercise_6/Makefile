test :
	- rm -rf *.out
	(time $(TEST_KELVIN) kelvin.conf) >kelvin.out 2>&1
	perl $(KELVIN_ROOT)/utils/wordDiff.pl br.out br.out-baseline 0
	diff -i ppl.out ppl.out-baseline

clean :
	- rm -f *.out *.so *memory.dat compiled
