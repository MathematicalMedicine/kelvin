# Originally checked against iterative output w/gnuplot of bayes ratio vs position

test :
	- rm -rf *.out
	(time $(TEST_KELVIN) kelvin.conf) >kelvin.out 2>&1
	diff -i br.out br.out-baseline

test-USE_DL :
	- rm -rf *.out *.so compiled
	(time $(KELVIN_ROOT)/kelvin-POLYCOMP_DL kelvin.conf) >kelvin-COMP_DL.out 2>&1
	diff -i br.out br.out-baseline
	(time $(KELVIN_ROOT)/kelvin-POLYUSE_DL kelvin.conf) >kelvin-USE_DL.out 2>&1
	diff -i br.out br.out-baseline

test-FIXED :
	-rm -f LRTest*
	(perl $(KELVIN_ROOT)/LRTest.pl) >kelvin-FIXED.out 2>&1

clean :
	- rm -rf *.out *.so *memory.dat compiled LRTest*

