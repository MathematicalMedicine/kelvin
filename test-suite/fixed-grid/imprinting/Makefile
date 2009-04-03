test :
	- rm -rf *.out
	(time $(TEST_KELVIN) kelvin.conf) >kelvin.out 2>&1
	diff -i br.out br.out-baseline
	diff -i ppl.out ppl.out-baseline

test-USE_DL :

clean :
	- rm -rf *.out *.so *memory.dat compiled
