# Identical to the SA_DT under dynamic PE

test-LKS :
	- rm -rf *.out adhoc* *memory.dat
	$(KELVIN_ROOT)/LKS/remove_study.sh
	(time $(TEST_KELVIN) kelvin.conf) >kelvin.log 2>&1 || log "Target $@ FAILED"
	(time ./cycle.sh) >cycle.out 2>&1 || log "Target $@ FAILED"
	perl $(KELVIN_ROOT)/utils/wordDiff.pl br.out br.out-baseline .000001

test-USE_DL :
	echo "Not available"

test-FIXED :
	echo "Not available"

clean :
	- rm -rf *.out adhoc* *memory.dat
