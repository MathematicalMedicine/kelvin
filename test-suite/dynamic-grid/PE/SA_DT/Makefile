# Originally checked against iterative output w/gnuplot of bayes ratio vs position

test :
	- rm -rf *.out
	(time $(TEST_KELVIN) kelvin.conf) >kelvin.out 2>&1 || log "$(TEST_KELVIN) run for target $@ in $(PWD) FAILED"
	perl $(KELVIN_ROOT)/utils/wordDiff.pl br.out br.out-baseline .000001 || log "wordDiff for target $@ in $(PWD) FAILED"

test-USE_DL :
	- rm -rf *.out *.so compiled *.linking *.compiling
	(time $(TEST_KELVIN)-POLYCOMP_DL kelvin.conf) >kelvin-COMP_DL.out 2>&1 || log "$(TEST_KELVIN) run for target $@ in $(PWD) FAILED"
	perl $(KELVIN_ROOT)/utils/wordDiff.pl br.out br.out-baseline .000001 || log "wordDiff for target $@ in $(PWD) FAILED"
	(time $(TEST_KELVIN)-POLYUSE_DL kelvin.conf) >kelvin-USE_DL.out 2>&1 || log "$(TEST_KELVIN) run for target $@ in $(PWD) FAILED"
	perl $(KELVIN_ROOT)/utils/wordDiff.pl br.out br.out-baseline .000001 || log "wordDiff for target $@ in $(PWD) FAILED"

test-FIXED :
	-rm -f LRTest*
	(perl $(KELVIN_ROOT)/LRTest.pl) >kelvin-FIXED.out 2>&1 || log "$(TEST_KELVIN) and $(KELVIN_ROOT)/LRTest.pl FIXED run for target $@ in $(PWD) FAILED"

clean :
	- rm -rf *.out *.so *memory.dat compiled *.linking *.compiling LRTest*

