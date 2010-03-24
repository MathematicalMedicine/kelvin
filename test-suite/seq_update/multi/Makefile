test :
	- rm -rf br.updated
	@echo "NOTE: This test is expected to warn about missing positions."
	@echo "      This warning is not a failure."
	$(SEQUPDATE_BINARY) -m br1.out br2.out br3.out > br.updated
	perl $(KELVIN_ROOT)/utils/wordDiff.pl br.updated br.updated-baseline .000001

clean :
	- rm -rf br.updated
