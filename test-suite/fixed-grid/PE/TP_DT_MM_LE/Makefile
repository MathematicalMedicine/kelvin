test :
	- rm -rf *.out
	(time $(TEST_KELVIN) kelvin.conf) >kelvin.log 2>&1 || log "Target $@ FAILED"
	perl $(KELVIN_ROOT)/utils/wordDiff.pl ppl.out ppl.out-baseline .000001 || log "Target $@ FAILED"

test-USE_DL :
	- rm -rf *.out *.so compiled

clean :
	- rm -rf *.out *.so *memory.dat compiled
