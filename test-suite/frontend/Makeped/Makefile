test :
	- rm -f *.out
	(time $(KELVIN_SCRIPT) kelvin.conf) >kelvin.out 2>&1
	perl $(KELVIN_ROOT)/utils/wordDiff.pl br.out br.out-baseline .000001 || log "wordDiff for target $@ in $(PWD) FAILED"
	perl $(KELVIN_ROOT)/utils/wordDiff.pl ppl.out ppl.out-baseline .000001 || log "wordDiff for target $@ in $(PWD) FAILED"

clean :
	- rm -f *.out

