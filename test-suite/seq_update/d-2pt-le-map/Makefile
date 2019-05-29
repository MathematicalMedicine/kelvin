test :
	- rm -rf br.updated ppl.updated
	$(SEQUPDATE_BINARY) -M map -O br.updated br1.out br2.out br3.out > ppl.updated
	perl $(KELVIN_ROOT)/utils/wordDiff.pl br.updated br.updated-baseline .000001 || log "wordDiff for target $@ in $(PWD) FAILED"
	perl $(KELVIN_ROOT)/utils/wordDiff.pl ppl.updated ppl.updated-baseline .000001 || log "wordDiff for target $@ in $(PWD) FAILED"

clean :
	- rm -rf br.updated ppl.updated
