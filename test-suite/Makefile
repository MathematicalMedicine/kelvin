# Makefile to perform kelvin tests

# find -maxdepth DNE on Solaris...
#TESTS = $(shell find * -maxdepth 0 -type d)
SHELL = bash
TESTS = $(shell ls -d */ | cut -f1 -d\/)

.PHONY : test clean $(TESTS)

check : ACTION := check
check : $(TESTS)

test : ACTION := test
test : $(TESTS)

test-USE_DL : ACTION := test-USE_DL
test-USE_DL : $(TESTS)

test-LKS : ACTION := test-LKS
test-LKS : $(TESTS)

test-FIXED : ACTION := test-FIXED
test-FIXED : $(TESTS)

clean : ACTION := clean
clean : $(TESTS)

$(TESTS) : 
	$(MAKE) -C $@ -f Makefile $(ACTION)
