test :
	- rm -f kelvin.freq *.out
	(time $(KELVIN_SCRIPT) kelvin.conf) >kelvin.log 2>&1
	diff -i kelvin.freq freq.dat-baseline

clean :
	- rm -f kelvin.freq *.out
