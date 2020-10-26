# Makefile to perform kelvin tests

# find -maxdepth DNE on Solaris...
#TESTS = $(shell find * -maxdepth 0 -type d)
SHELL = bash
TESTS = $(shell ls -d */ | cut -f1 -d\/)

.PHONY : test test-USE_DL clean $(TESTS)

test : ACTION := test
test : $(TESTS)

test-USE_DL :

check :

test-FIXED :

clean :

$(TESTS) : 
#	($(KELVIN_ROOT)/kdiff.pl -p1 $@ -p2 $@ -c1 /dev/null -c2 /dev/null)
	rm $@/kdiff.out 2>/dev/null || true
	($(KELVIN_ROOT)/kdiff.pl -p1 original -p2 $@ -c1 /dev/null -c2 /dev/null) >$@/kdiff.out 2>&1 || true
	diff $@/kdiff.out $@/kdiff.out-baseline || log "Target $@ FAILED"
