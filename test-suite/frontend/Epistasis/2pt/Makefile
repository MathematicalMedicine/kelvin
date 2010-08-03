test :
	- rm -rf *.out
	(time $(KELVIN_SCRIPT) kelvin.conf) >kelvin.out 2>&1
	perl $(KELVIN_ROOT)/utils/wordDiff.pl br.out br.out-baseline 0
	perl $(KELVIN_ROOT)/utils/wordDiff.pl ppl.out ppl.out-baseline 0
	perl $(KELVIN_ROOT)/utils/wordDiff.pl mod.out mod.out-baseline 0

clean :
	- rm -rf *.out *.pass?
