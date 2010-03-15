TEST_UPDATE=../../../seq_update/calc_updated_ppl

test :
	- rm -rf br.updated ppl.updated
	$(TEST_UPDATE) -M map -O br.updated br1.out br2.out br3.out > ppl.updated
	perl $(KELVIN_ROOT)/utils/wordDiff.pl br.updated br.updated-baseline .000001
	perl $(KELVIN_ROOT)/utils/wordDiff.pl ppl.updated ppl.updated-baseline .000001

clean :
	- rm -rf br.updated ppl.updated
