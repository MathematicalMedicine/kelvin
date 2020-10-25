test :
	- rm -f kelvin.freq *.out *.tmp
	(time $(KELVIN_SCRIPT) kelvin.conf) >kelvin.log 2>&1
	diff -i kelvin.freq freq.dat-baseline
	perl $(KELVIN_ROOT)/utils/wordDiff.pl br.out br.out-baseline 0.000001 || log "wordDiff for target $@ in $(PWD) FAILED"
	perl $(KELVIN_ROOT)/utils/wordDiff.pl ppl.out ppl.out-baseline 0.000001 || log "wordDiff for target $@ in $(PWD) FAILED"

clean :
	- rm -f kelvin.freq *.out *.tmp
