SEQUPDATE_BINARY ?= missing

test :
	- rm -rf br.updated
	$(SEQUPDATE_BINARY) -P ppl.out -M map.dat -f -U br1.out br2.out > br.updated
	perl $(KELVIN_ROOT)/utils/wordDiff.pl br.updated br.updated-baseline .000001 || log "Target $@ FAILED"

clean :
	- rm -rf br.updated 
