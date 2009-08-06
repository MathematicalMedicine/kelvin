test :
	- rm -rf *.out
	(time $(TEST_KELVIN) kelvin.conf) >kelvin.out 2>&1
	diff -i br.out br.out-baseline
	diff -i ppl.out ppl.out-baseline

test-FIXED :
	-rm -f LRTest*
	perl $(KELVIN_ROOT)/LRTest.pl

clean :
	- rm -rf *.out *.so *memory.dat compiled
