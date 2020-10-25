test :
	- rm -f br.out ppl.out kelvin.log
	(time $(TEST_KELVIN) kelvin.conf) >kelvin.log 2>&1 || log "$(TEST_KELVIN) run for target $@ in $(PWD) FAILED"
	perl $(KELVIN_ROOT)/utils/wordDiff.pl br.out br.out-baseline .000001 || log "wordDiff for target $@ in $(PWD) FAILED"
	perl $(KELVIN_ROOT)/utils/wordDiff.pl ppl.out ppl.out-baseline .000001 || log "wordDiff for target $@ in $(PWD) FAILED"

test-FIXED :
	echo Fixed-grid comparison not available

test-USE_DL :
	echo Compilation not appropriate for this test.

clean :
	- rm -rf *.so *.out *memory.dat *~ compiled
