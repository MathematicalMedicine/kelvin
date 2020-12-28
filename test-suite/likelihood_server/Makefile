# Makefile to perform kelvin tests

# find -maxdepth DNE on Solaris...
#TESTS = $(shell find * -maxdepth 0 -type d)
SHELL = bash
TESTS = $(shell ls -d */ | cut -f1 -d\/)

.PHONY : test-LKS clean $(TESTS)

check : ACTION := check
check : $(TESTS)

test-LKS : ACTION := test-LKS
test-LKS : $(TESTS)

test-FIXED :

test-USE_DL :

clean : ACTION := clean
clean : $(TESTS)

$(TESTS) : 
	$(MAKE) -C $@ -f Makefile $(ACTION)

