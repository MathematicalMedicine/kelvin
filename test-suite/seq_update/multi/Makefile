TEST_UPDATE=../../../seq_update/calc_updated_ppl

test :
	- rm -rf br.updated
	@echo "NOTE: This test is expected to warn about missing positions."
	@echo "      This warning is not a failure."
	$(TEST_UPDATE) -m br1.out br2.out br3.out > br.updated
	diff -i br.updated br.updated-baseline

clean :
	- rm -rf br.updated
