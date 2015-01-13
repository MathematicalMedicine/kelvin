# Identical to the SA_DT under dynamic PE

test :
	- rm -rf *.out adhoc* *memory.dat
	$(KELVIN_ROOT)/LKS/remove_study.sh 4
	(time /export/local/bcmmtools/LKS/merlin_lk_cycle.sh) >cycle.out 2>&1
	perl $(KELVIN_ROOT)/utils/wordDiff.pl br.out br.out-baseline .000001

test-USE_DL :
	echo "Not available"

test-FIXED :
	echo "Not available"

clean :
	- rm -rf *.out adhoc* *memory.dat
