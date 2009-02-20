test :
	- rm -f br.out ppl.out kelvin.out
	time $(TEST_KELVIN) kelvin.conf >kelvin.out 2>&1
	perl $(KELVIN_ROOT)/utils/wordDiff.pl br.out br.out-baseline .000001
	diff -i ppl.out ppl.out-baseline

clean :
	- rm -f *.so *.out *memory.dat *~
