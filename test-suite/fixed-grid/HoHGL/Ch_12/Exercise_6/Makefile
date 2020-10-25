test :
	- rm -rf *.out
	(time $(TEST_KELVIN) kelvin.conf) >kelvin.log 2>&1 || log "$(TEST_KELVIN) run for target $@ in $(PWD) FAILED"
	perl $(KELVIN_ROOT)/utils/wordDiff.pl br.out br.out-baseline 0 || log "wordDiff for target $@ in $(PWD) FAILED"
	perl $(KELVIN_ROOT)/utils/wordDiff.pl ppl.out ppl.out-baseline .000001 || log "wordDiff for target $@ in $(PWD) FAILED"

clean :
	- rm -f *.out *.so *memory.dat compiled
