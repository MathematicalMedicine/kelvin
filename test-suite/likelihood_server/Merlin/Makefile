# Identical to the SA_DT under dynamic PE

test-LKS :
	- rm -rf *.out adhoc* *memory.dat
	$(KELVIN_ROOT)/LKS/remove_study.sh
	(time $(TEST_KELVIN) kelvin.conf) >kelvin.log 2>&1 || log "Target $@ FAILED"
	cut -d ' ' -f 1-5 br.out-baseline >br.out-baseline.first_5
	(time ./cycle.sh) >cycle.out 2>&1 || log "Target $@ FAILED"
	cut -d ' ' -f 1-5 br.out >br.out.first_5
	perl $(KELVIN_ROOT)/utils/wordDiff.pl br.out.first_5 br.out-baseline.first_5 .0 || log "Target $@ FAILED"

test-USE_DL :
	echo "Not available"

test-FIXED :
	echo "Not available"

clean :
	- rm -rf *.out adhoc* *memory.dat
