test :
	- rm -rf single-* dynamic-* fixed-* *.out
	# Produce the reference using dynamic grid
	$(TEST_KELVIN) kelvin.conf >kelvin.out 2>&1
	# Compare dynamic and fixed for each marker
	perl -n mod2singles.pl <mod.out

clean :
	- rm -rf single-* dynamic-* fixed-* *.out
