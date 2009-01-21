# Makefile to perform kelvin tests

# find -maxdepth DNE on Solaris...
#TESTS = $(shell find * -maxdepth 0 -type d)
SHELL = bash
TEST_KELVIN = ~/kelvin/trunk/kelvin
KELVIN_ROOT = ~/kelvin/trunk/
export TEST_KELVIN KELVIN_ROOT
TESTS = $(shell ls -d */ | cut -f1 -d\/)

.PHONY : test clean $(TESTS)

test : ACTION := test
test : $(TESTS)

test-USE_DL : ACTION := test-USE_DL
test-USE_DL : $(TESTS)

clean : ACTION := clean
clean : $(TESTS)

$(TESTS) : 
	$(MAKE) -C $@ -f Makefile $(ACTION)
