# Makefile to perform kelvin tests

TESTS = \
	consanguinity_loop \
	marriage_loop \
	broken_loop \
	multi-founder \
	imprinting \
        PE \
        NonPE \
	HoHGL

.PHONY : all $(TESTS)

all : $(TESTS)

$(TESTS) : 
	$(MAKE) -C $@ -f Makefile test

