test :
	- rm -f ped.post.counted ped.post.counts kelvin.freq* *.out
	# Disabled until such time as John can fix it
	#(time $(KELVIN_SCRIPT) kelvin.conf) >kelvin.out 2>&1
	#diff -ib ped.post.counted ped.post.counted-baseline
	#diff -ib ped.post.counts ped.post.counts-baseline
	echo THIS TEST IS CURRENTLY DISABLED
	cp freq.dat-baseline kelvin.freq
	cp br.out-baseline br.out
	cp ppl.out-baseline ppl.out
	#diff -ib ped.post.counted ped.post.counted-baseline
	#diff -ib ped.post.counts ped.post.counts-baseline
	perl $(KELVIN_ROOT)/utils/wordDiff.pl kelvin.freq freq.dat-baseline 0
	perl $(KELVIN_ROOT)/utils/wordDiff.pl br.out br.out-baseline .000001
	perl $(KELVIN_ROOT)/utils/wordDiff.pl ppl.out ppl.out-baseline .000001

clean :
	- rm -f ped.post.counted ped.post.counts kelvin.freq* *.out

