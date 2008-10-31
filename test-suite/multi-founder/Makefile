test :
	- rm -f br.out ppl.out kelvin.out
	time ../../kelvin kelvin.conf >kelvin.out 2>&1
	diff -i br.out br.out-baseline
	diff -i ppl.out ppl.out-baseline

clean :
	- rm -f *.so *.out *memory.dat *~
