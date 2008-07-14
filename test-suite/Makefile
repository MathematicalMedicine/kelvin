# Makefile to perform kelvin tests

TESTS = \
	cosanguinous_loop \
	marriage_loop \
        PE \
        NonPE

.PHONY : all $(TESTS)

all : $(TESTS)

$(TESTS) : 
	$(MAKE) -C $@ -f Makefile test

