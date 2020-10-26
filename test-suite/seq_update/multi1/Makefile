SEQUPDATE_BINARY ?= missing

test :
	- rm -rf br.updated
	$(SEQUPDATE_BINARY) -m br1.out br2.out br3.out > br.updated
	perl $(KELVIN_ROOT)/utils/wordDiff.pl br.updated br.updated-baseline .000001 || log "Target $@ FAILED"

clean :
	- rm -rf br.updated
