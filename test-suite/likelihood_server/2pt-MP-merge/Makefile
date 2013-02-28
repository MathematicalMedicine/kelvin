# Identical to the SA_DT under dynamic PE

test :
	- rm -rf *.out adhoc* *memory.dat study_*.dat
	$(KELVIN_ROOT)/LKS/remove_study.sh
	(time ./cycle.sh) >cycle.out 2>&1
	perl $(KELVIN_ROOT)/utils/wordDiff.pl br.out-all-pedigrees br.out-baseline .000001

test-USE_DL :
	echo "Not available"

test-FIXED :
	echo "Not available"

clean :
	- rm -rf *.out adhoc* *memory.dat study_*.dat
