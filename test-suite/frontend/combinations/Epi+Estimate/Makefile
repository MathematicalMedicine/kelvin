test :
	- rm -f *.out kelvin.freq *pass?
	(time $(KELVIN_SCRIPT) kelvin.conf) >kelvin.log 2>&1
	diff -ib kelvin.freq-baseline kelvin.freq
	diff -ib ped.post.pass0-baseline ped.post.pass0
	diff -ib ped.post.pass1-baseline ped.post.pass1
	diff -ib locus.dat.pass0-baseline locus.dat.pass0
	diff -ib locus.dat.pass1-baseline locus.dat.pass1

clean :
	- rm -f *.out kelvin.freq *pass?

