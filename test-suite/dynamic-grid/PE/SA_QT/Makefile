# Originally checked against iterative output w/gnuplot of bayes ratio vs position

test :
	- rm -f *.out
	(time $(TEST_KELVIN) kelvin.conf) >kelvin.log 2>&1 || log "Target $@ FAILED"
	perl $(KELVIN_ROOT)/utils/wordDiff.pl br.out br.out-baseline .000001 || log "Target $@ FAILED"

#test-FIXED :
#	-rm -f LRTest*
#	(perl $(KELVIN_ROOT)/LRTest.pl) >kelvin-FIXED.log 2>&1
test-FIXED :
	- rm -rf *.out
	@echo -e "\n\
    **** THIS TEST IS CURRENTLY DISABLED AS QTMeanMode and QTStandardDevMode are not compatible with FixedModels ****\n\
"


test-USE_DL :
	- rm -rf *.out *.so compiled *.linking *.compiling
	(time $(TEST_KELVIN)-POLYCOMP_DL kelvin.conf) >kelvin-COMP_DL.log 2>&1 || log "Target $@ FAILED"
	perl $(KELVIN_ROOT)/utils/wordDiff.pl br.out br.out-baseline .000001 || log "Target $@ FAILED"
	(time $(TEST_KELVIN)-POLYUSE_DL kelvin.conf) >kelvin-USE_DL.log 2>&1 || log "Target $@ FAILED"
	perl $(KELVIN_ROOT)/utils/wordDiff.pl br.out br.out-baseline .000001 || log "Target $@ FAILED"

clean :
	- rm -rf *.out *.so *memory.dat compiled *.linking *.compiling LRTest*

