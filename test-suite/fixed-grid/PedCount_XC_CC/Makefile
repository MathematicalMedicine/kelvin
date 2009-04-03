test :
	- rm -rf PC* *.out
	$(KELVIN_ROOT)/PedCount.pl -config kelvin.conf -count -XC -write
	(time $(TEST_KELVIN) kelvin-counts.conf) >kelvin.out 2>&1
	diff -i br.out br.out-baseline
	diff -i ppl.out ppl.out-baseline

clean :
	- rm -rf PC* *.out
