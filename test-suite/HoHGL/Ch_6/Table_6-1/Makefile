test :
	- rm -f br.out ppl.out kelvin.out
	time ../../../../kelvin kelvin.conf >kelvin.out 2>&1
	perl ../../../../utils/wordDiff.pl br.out br.out-baseline 0
	diff -i ppl.out ppl.out-baseline

clean :
	- rm -f *.out *memory.dat *~
