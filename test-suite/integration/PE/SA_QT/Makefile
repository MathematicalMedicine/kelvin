# Originally checked against iterative output w/gnuplot of bayes ratio vs position

test :
	- rm -f br.out ppl.out kelvin.out
	time $(TEST_KELVIN) kelvin.conf >kelvin.out 2>&1
	perl $(KELVIN_ROOT)/utils/wordDiff.pl br.out br.out-baseline .000001

test-USE_DL :
	echo Compilation not available for QT runs yet
#	- rm -rf br.out ppl.out kelvin-*_DL.out *.so compiled
#	time $(KELVIN_ROOT)/kelvin-POLYCOMP_DL kelvin.conf >kelvin-COMP_DL.out 2>&1
#	perl $(KELVIN_ROOT)/utils/wordDiff.pl br.out br.out-baseline .000001
#	time $(KELVIN_ROOT)/kelvin-POLYUSE_DL kelvin.conf >kelvin-USE_DL.out 2>&1
#	perl $(KELVIN_ROOT)/utils/wordDiff.pl br.out br.out-baseline .000001

clean :
	- rm -f *.so *.out *memory.dat *~

