test :
	- rm -f br.out ppl.out kelvin.out
	- time ../../kelvin kelvin.conf >&kelvin.out
	grep "Not all loops in pedigrees are broken" kelvin.out
