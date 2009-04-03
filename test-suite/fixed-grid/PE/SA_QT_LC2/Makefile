test :
	- rm -rf *.out
	(time $(TEST_KELVIN) kelvin.conf) >kelvin.out 2>&1
	perl $(KELVIN_ROOT)/utils/wordDiff.pl br.out br.out-baseline .000001

test-USE_DL :
	echo Compilation not available for QT runs yet

clean :
	- rm -rf *.out *.so *memory.dat compiled

