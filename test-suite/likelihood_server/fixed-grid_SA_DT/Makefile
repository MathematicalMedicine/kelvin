# Identical to the SA_DT under dynamic PE

test-LKS :
	- rm -rf *.out adhoc* *memory.dat study_*.dat
	$(KELVIN_ROOT)/LKS/remove_study.sh
	(time $(KELVIN_ROOT)/LKS/cycle.sh) >cycle.out 2>&1 || log "Target $@ FAILED"
	perl $(KELVIN_ROOT)/utils/wordDiff.pl br.out br.out-baseline .000001 || log "Target $@ FAILED"

test-USE_DL :
	echo "Not available"

test-FIXED :
	echo "Not available"

clean :
	- rm -rf *.out adhoc* *memory.dat study_*.dat
