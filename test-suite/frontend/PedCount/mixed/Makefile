test :
	- rm -rf PC* *.out
	$(KELVIN_ROOT)/PedCount.pl -config kelvin.conf -count -write
# Write a new baseline when needed, but this is slow
#	(time $(TEST_KELVIN) kelvin.conf) >kelvin.out 2>&1
	(time $(TEST_KELVIN) kelvin-counts.conf) >kelvin-counts.out 2>&1
	perl $(KELVIN_ROOT)/utils/wordDiff.pl br.out br.out-baseline .000001
	perl $(KELVIN_ROOT)/utils/wordDiff.pl ppl.out ppl.out-baseline .000001

test-USE_DL :

clean :
	- rm -rf PC* *.out
