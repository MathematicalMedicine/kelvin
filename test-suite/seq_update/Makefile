# Makefile to perform kelvin tests

# find -maxdepth DNE on Solaris...
#TESTS = $(shell find * -maxdepth 0 -type d)
SHELL = bash
TESTS = $(shell ls -d */ | cut -f1 -d\/)

.PHONY : test test-USE_DL clean $(TESTS)

test : ACTION := test
test : $(TESTS)

check :
	$(MAKE) -C d-2pt-le -f Makefile test

test-USE_DL : 

test-FIXED :

clean : ACTION := clean
clean : $(TESTS)

$(TESTS) : 
	$(MAKE) -C $@ -f Makefile $(ACTION)
