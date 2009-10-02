TEST_UPDATE=../../../seq_update/calc_updated_ppl

test :
	- rm -rf br.updated
	$(TEST_UPDATE) -m br1.out br2.out br3.out > br.updated
	diff -i br.updated br.updated-baseline

clean :
	- rm -rf br.updated
