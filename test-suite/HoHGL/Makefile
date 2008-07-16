# Makefile to perform kelvin tests

TESTS = \
	Exercises_3_and_5 \
        Exercise_6

.PHONY : test $(TESTS)

test : $(TESTS)

$(TESTS) : 
	$(MAKE) -C $@ -f Makefile test

