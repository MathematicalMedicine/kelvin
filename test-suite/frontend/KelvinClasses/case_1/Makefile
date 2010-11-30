test :
	- rm -rf *.new
	perl -I $(KELVIN_ROOT) classtest.pl
	diff -b ped.new ped.post-test
	diff -b locus.new locus.dat-test
	diff -b map.new map.dat-test
	diff -b freq.new freq.dat-test

clean :
	- rm -rf *.new
