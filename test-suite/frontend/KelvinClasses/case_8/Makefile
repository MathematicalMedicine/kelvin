test :
	- rm -rf report.txt
	perl -I $(KELVIN_ROOT) classtest.pl
	diff -bB report.txt report.txt-baseline
	- rm -rf report.txt

clean :
	- rm -rf report.txt
