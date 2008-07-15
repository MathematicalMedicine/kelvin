# Makefile to perform kelvin tests

TESTS = \
	cosanguinous_loop \
	marriage_loop \
	broken_loop \
	multi-founder \
        PE \
        NonPE

.PHONY : all $(TESTS)

all : $(TESTS)

$(TESTS) : 
	$(MAKE) -C $@ -f Makefile test

