# Makefile to perform kelvin tests

TESTS = $(shell find * -maxdepth 0 -type d)

.PHONY : test clean $(TESTS)

test : ACTION := test
test : $(TESTS)

clean : ACTION := clean
clean : $(TESTS)

$(TESTS) : 
	$(MAKE) -C $@ -f Makefile $(ACTION)
