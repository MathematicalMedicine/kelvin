test :
	- rm -f kelvin.freq *.out
	(time $(KELVIN_SCRIPT) kelvin.conf --SkipPedCount) >kelvin.out 2>&1
	diff -i kelvin.freq freq.dat-baseline
	perl $(KELVIN_ROOT)/utils wordDiff.pl br.out br.out-baseline 0.000001
	perl $(KELVIN_ROOT)/utils wordDiff.pl ppl.out ppl.out-baseline 0.000001
	(time $(KELVIN_SCRIPT) kelvin.conf) >kelvin.out 2>&1 
	perl $(KELVIN_ROOT)/utils wordDiff.pl br.out br.out-baseline 0.000001
	perl $(KELVIN_ROOT)/utils wordDiff.pl ppl.out ppl.out-baseline 0.000001

clean :
	- rm -f kelvin.freq *.out
