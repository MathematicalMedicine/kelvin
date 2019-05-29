test :
	- rm -rf *.out
	(time $(KELVIN_SCRIPT) kelvin.conf) >kelvin.out 2>&1
	perl $(KELVIN_ROOT)/utils/wordDiff.pl br.out br.out-baseline 0 || log "wordDiff for target $@ in $(PWD) FAILED"
	perl $(KELVIN_ROOT)/utils/wordDiff.pl mod.out mod.out-baseline 0 || log "wordDiff for target $@ in $(PWD) FAILED"

clean :
	- rm -rf *.out *.pass?
