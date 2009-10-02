TEST_UPDATE=../../../seq_update/calc_updated_ppl

test :
	- rm -rf br.updated ppl.updated
	$(TEST_UPDATE) -M map -O br.updated br1.out br2.out br3.out > ppl.updated
	diff -i br.updated br.updated-baseline
	diff -i ppl.updated ppl.updated-baseline

clean :
	- rm -rf br.updated ppl.updated
