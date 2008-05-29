test :
	- rm -f avghet.out ppl.out kelvin.out
	time ../../kelvin kelvin.conf >&kelvin.out
	diff avghet.out avghet.out-baseline
	diff ppl.out ppl.out-baseline
