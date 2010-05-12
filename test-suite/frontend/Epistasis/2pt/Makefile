test :
	- rm -rf *.out
	(time $(KELVIN_SCRIPT) kelvin.conf) >kelvin.out 2>&1
	diff -I '^# Version' -i br.out br.out-baseline
	diff -I '^# Version' -i ppl.out ppl.out-baseline
	diff -I '^# Version' -i mod.out mod.out-baseline

clean :
	- rm -rf *.out *.so *memory.dat compiled
