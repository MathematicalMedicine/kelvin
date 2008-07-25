test :
	- rm -f br.out ppl.out kelvin.out
	time ../../../kelvin kelvin.conf >&kelvin.out
	- diff br.out br.out-baseline
	diff ppl.out ppl.out-baseline
