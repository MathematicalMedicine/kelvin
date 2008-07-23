# Makefile to perform kelvin tests

TESTS = \
	Ch_6 \
        Ch_12

.PHONY : test $(TESTS)

test : $(TESTS)

$(TESTS) : 
	$(MAKE) -C $@ -f Makefile test

