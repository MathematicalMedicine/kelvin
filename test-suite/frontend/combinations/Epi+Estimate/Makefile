test :
	- rm -f *.out kelvin.freq* *pass?
	(time $(KELVIN_SCRIPT) kelvin.conf) >kelvin.out 2>&1
	#diff -ib ped.post.counted ped.post.counted-baseline
	#diff -ib ped.post.counts ped.post.counts-baseline
	diff -ib freq.dat-baseline kelvin.freq
	perl $(KELVIN_ROOT)/utils/wordDiff.pl br.out br.out-baseline .000001
	perl $(KELVIN_ROOT)/utils/wordDiff.pl ppl.out ppl.out-baseline .000001

clean :
	- rm -f *.out kelvin.freq* *pass?

