# Makefile to perform kelvin tests

TESTS = \
	Table_6-1 \
        Table_6-2

.PHONY : test $(TESTS)

test : $(TESTS)

$(TESTS) : 
	$(MAKE) -C $@ -f Makefile test

