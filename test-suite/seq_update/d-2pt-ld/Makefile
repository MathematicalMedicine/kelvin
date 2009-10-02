TEST_UPDATE=../../../seq_update/calc_updated_ppl

test :
	- rm -rf br.updated ppl.updated
	$(TEST_UPDATE) -O br.updated br1.out br2.out > ppl.updated
	diff -ib br.updated br.updated-baseline
	diff -ib ppl.updated ppl.updated-baseline

clean :
	- rm -rf br.updated ppl.updated
