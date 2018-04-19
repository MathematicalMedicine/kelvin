test :
	- rm -rf *.new
	perl -I $(KELVIN_ROOT) classtest-1.pl
	diff -b locus.new locus.dat-baseline
	diff -b freq.new freq.dat-baseline
	perl -I $(KELVIN_ROOT) classtest-2.pl
	diff -b freq.new freq.dat-baseline
	diff -b map.new map.dat-baseline
	perl -I $(KELVIN_ROOT) classtest-3.pl
	diff -b locus.new locus.dat-baseline
	diff -b map.new map.dat-baseline
	perl -I $(KELVIN_ROOT) classtest-4.pl
	diff -b freq.new freq.dat-baseline
	diff -b map.new map.dat-baseline
	diff -b locus.new locus.dat-baseline
	- rm -rf *.new

clean :
	- rm -rf *.new
