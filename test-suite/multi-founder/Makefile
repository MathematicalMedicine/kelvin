test :
	- rm -f br.out ppl.out kelvin.out
	time ../../kelvin kelvin.conf >kelvin.out 2>&1
	diff br.out br.out-baseline
	diff ppl.out ppl.out-baseline

clean :
	- rm -f *.out *memory.dat *~
