test :
	- rm -rf *.out
	(time $(TEST_KELVIN) kelvin.conf) >kelvin.out 2>&1 || log "$(TEST_KELVIN) run for target $@ in $(PWD) FAILED"
	perl $(KELVIN_ROOT)/utils/wordDiff.pl br.out br.out-baseline .000001 || log "wordDiff for target $@ in $(PWD) FAILED"

clean :
	- rm -rf *.out *.so *memory.dat compiled
