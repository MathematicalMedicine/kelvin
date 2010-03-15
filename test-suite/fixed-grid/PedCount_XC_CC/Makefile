test :
	- rm -rf PC* *.out
	$(KELVIN_ROOT)/PedCount.pl -config kelvin.conf -count -XC -write
	(time $(TEST_KELVIN) kelvin-counts.conf) >kelvin.out 2>&1
	perl $(KELVIN_ROOT)/utils/wordDiff.pl br.out br.out-baseline .000001
	perl $(KELVIN_ROOT)/utils/wordDiff.pl ppl.out ppl.out-baseline .000001

test-USE_DL :

clean :
	- rm -rf PC* *.out
