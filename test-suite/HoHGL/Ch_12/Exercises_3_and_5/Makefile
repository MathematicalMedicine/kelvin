test :
	- rm -f br.out ppl.out kelvin.out
	time ../../../../kelvin kelvin.EX3 >kelvin.out-EX3 2>&1
	diff -i br.out br.out-EX3
	diff -i ppl.out ppl.out-EX3
	time ../../../../kelvin kelvin.EX5a >kelvin.out-EX5a 2>&1
	diff -i br.out br.out-EX5a
	diff -i ppl.out ppl.out-EX5a
	time ../../../../kelvin kelvin.EX5b >kelvin.out-EX5b 2>&1
	diff -i br.out br.out-EX5b
	diff -i ppl.out ppl.out-EX5b
	time ../../../../kelvin kelvin.EX5c >kelvin.out-EX5c 2>&1
	diff -i br.out br.out-EX5c
	diff -i ppl.out ppl.out-EX5c

clean :
	- rm -f *.so *.out kelvin.out-EX* *memory.dat *~

