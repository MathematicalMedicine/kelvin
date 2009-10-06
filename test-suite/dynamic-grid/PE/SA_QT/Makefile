# Originally checked against iterative output w/gnuplot of bayes ratio vs position

test :
	- rm -f *.out
	(time $(TEST_KELVIN) kelvin.conf) >kelvin.out 2>&1
	perl $(KELVIN_ROOT)/utils/wordDiff.pl br.out br.out-baseline .000001

test-FIXED :
	-rm -f LRTest*
	(perl $(KELVIN_ROOT)/LRTest.pl) >kelvin-FIXED.out 2>&1

clean :
	- rm -rf *.out *.so *memory.dat compiled LRTest*

