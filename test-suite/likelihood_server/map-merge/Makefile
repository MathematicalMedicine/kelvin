test-LKS :
	- rm -rf *.out adhoc* *memory.dat study_*.dat
	$(KELVIN_ROOT)/LKS/remove_study.sh
	(time ./cycle.sh) >cycle.out 2>&1
	perl $(KELVIN_ROOT)/utils/wordDiff.pl br.out br.out-baseline .000001 || log "Target $@ FAILED"

test-USE_DL :
	echo "Not available"

test-FIXED :
	echo "Not available"

clean :
	- rm -rf *.out adhoc* *memory.dat study_*.dat
