# Makefile to perform kelvin tests

# find -maxdepth DNE on Solaris...
#TESTS = $(shell find * -maxdepth 0 -type d)
SHELL = bash
TESTS = $(shell ls -d */ | cut -f1 -d\/)

.PHONY : test clean $(TESTS)

test : ACTION := test
test : $(TESTS)

clean : ACTION := clean
clean : $(TESTS)

$(TESTS) : 
	$(MAKE) -C $@ -f Makefile $(ACTION)
