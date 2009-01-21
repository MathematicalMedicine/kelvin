test :
	- rm -f br.out ppl.out kelvin.out
	time $(TEST_KELVIN) kelvin.conf >kelvin.out 2>&1
	diff -i br.out br.out-baseline

test-USE_DL :
	- rm -rf br.out ppl.out kelvin-*_DL.out *.so compiled
	time ../../../kelvin-POLYCOMP_DL kelvin.conf >kelvin-COMP_DL.out 2>&1
	diff -i br.out br.out-baseline
	time ../../../kelvin-POLYUSE_DL kelvin.conf >kelvin-USE_DL.out 2>&1
	diff -i br.out br.out-baseline

clean :
	- rm -rf *.so *.out *memory.dat *~ compiled

