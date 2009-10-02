test :
	- rm -f br.out ppl.out kelvin.out
	(time $(TEST_KELVIN) kelvin.conf) >kelvin.out 2>&1
	diff -i br.out br.out-baseline
	diff -i ppl.out ppl.out-baseline

test-FIXED :
	echo Fixed-grid comparison not available

test-USE_DL :
	echo Compilation not appropriate for this test.

clean :
	- rm -rf *.so *.out *memory.dat *~ compiled
