test :
	- rm -rf *.out
	(time $(TEST_KELVIN) kelvin.conf) >kelvin.log 2>&1 || log "Target $@ FAILED"
	perl $(KELVIN_ROOT)/utils/wordDiff.pl br.out br.out-baseline .00001 || log "Target $@ FAILED"

test-USE_DL :
	- rm -rf *.out *.so compiled
	(time $(TEST_KELVIN)-POLYCOMP_DL kelvin.conf) >kelvin-COMP_DL.log 2>&1 || log "Target $@ FAILED"
	perl $(KELVIN_ROOT)/utils/wordDiff.pl br.out br.out-baseline .000001 || log "Target $@ FAILED"
	(time $(TEST_KELVIN)-POLYUSE_DL kelvin.conf) >kelvin-USE_DL.log 2>&1 || log "Target $@ FAILED"
	perl $(KELVIN_ROOT)/utils/wordDiff.pl br.out br.out-baseline .000001 || log "Target $@ FAILED"

clean :
	- rm -rf *.out *.so *memory.dat compiled
