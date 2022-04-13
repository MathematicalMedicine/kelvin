test :
	- rm -rf *.out
	(time $(KELVIN_SCRIPT) kelvin.conf) >kelvin.log 2>&1
	perl $(KELVIN_ROOT)/utils/wordDiff.pl br.out br.out-baseline 0.0001 || log "Target $@ FAILED"
	perl $(KELVIN_ROOT)/utils/wordDiff.pl mod.out mod.out-baseline 0.0001 || log "Target $@ FAILED"

clean :
	- rm -rf *.out *.pass?
