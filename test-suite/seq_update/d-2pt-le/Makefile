SEQUPDATE_BINARY ?= missing

test :
	- rm -rf br.updated ppl.updated
	$(SEQUPDATE_BINARY) -O br.updated br1.out br2.out br3.out > ppl.updated
	perl $(KELVIN_ROOT)/utils/wordDiff.pl br.updated br.updated-baseline .000001 || log "Target $@ FAILED"
	perl $(KELVIN_ROOT)/utils/wordDiff.pl ppl.updated ppl.updated-baseline .000001 || log "Target $@ FAILED"

clean :
	- rm -rf br.updated ppl.updated
