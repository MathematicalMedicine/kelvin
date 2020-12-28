# Makefile to perform kelvin tests

# find -maxdepth DNE on Solaris...
#TESTS = $(shell find * -maxdepth 0 -type d)
SHELL = bash
TESTS = $(shell ls -d */ | cut -f1 -d\/)

.PHONY : test test-USE_DL clean $(TESTS)

#test : ACTION := test
#test : $(TESTS)
test : 
	echo These tests have been disabled

check :
	$(MAKE) -C No_config -f Makefile test

test-USE_DL : ACTION := test-USE_DL
test-USE_DL : $(TESTS)

test-LKS :

test-FIXED :

clean : ACTION := clean
clean : $(TESTS)

$(TESTS) : 
	$(MAKE) -C $@ -f Makefile $(ACTION)
