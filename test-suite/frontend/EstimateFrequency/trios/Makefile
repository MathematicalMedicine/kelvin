test :
	- rm -f freq.dat *.out
	(time $(KELVIN_SCRIPT) kelvin.conf) >kelvin.out 2>&1
	diff -i freq.dat freq.dat-baseline

clean :
	- rm -f freq.dat *.out
