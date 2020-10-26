test :
	- rm -f kelvin.freq* *.out
	(time $(KELVIN_SCRIPT) kelvin.conf --SkipAnalysis) >kelvin.log 2>&1
	diff -i kelvin.freq freq.dat-baseline
	(time $(KELVIN_SCRIPT) kelvin.conf) >kelvin.log 2>&1 
	perl $(KELVIN_ROOT)/utils/wordDiff.pl br.out br.out-baseline 0.000001 || log "Target $@ FAILED"
	perl $(KELVIN_ROOT)/utils/wordDiff.pl ppl.out ppl.out-baseline 0.000001 || log "Target $@ FAILED"

clean :
	- rm -f kelvin.freq* *.out
