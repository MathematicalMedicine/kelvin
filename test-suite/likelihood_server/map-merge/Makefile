test :
	- rm -rf *.out adhoc* *memory.dat
	~/bcmmtools/remove_study.sh
	(time ./cycle.sh) >cycle.out 2>&1
	perl $(KELVIN_ROOT)/utils/wordDiff.pl br.out br.out-baseline .000001

test-USE_DL :
	echo "Not available"

test-FIXED :
	echo "Not available"

clean :
	- rm -rf *.out adhoc* *memory.dat
