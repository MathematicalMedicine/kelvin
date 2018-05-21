VERSION := $(shell echo `cat $(KELVIN_ROOT)/.maj`.`cat $(KELVIN_ROOT)/.min`.`cat $(KELVIN_ROOT)/.pat`)
test :
	- rm -rf *.out *.so compiled *.linking *.compiling
	(time $(KELVIN_ROOT)/kelvin-$(VERSION)-POLYCOMP_DL kelvin.conf) >kelvin-COMP_DL.out 2>&1
	perl $(KELVIN_ROOT)/utils/wordDiff.pl br.out br.out-baseline .000001
	(time $(KELVIN_ROOT)/kelvin-$(VERSION)-POLYUSE_DL kelvin.conf) >kelvin-USE_DL.out 2>&1
	perl $(KELVIN_ROOT)/utils/wordDiff.pl br.out br.out-baseline .000001

clean :
	- rm -rf *.out *.so *memory.dat compiled *.linking *.compiling LRTest*
