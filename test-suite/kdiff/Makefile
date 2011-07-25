test :
	$(KELVIN_ROOT)/kdiff.pl --ped --c1 /dev/null --c2 /dev/null --path1 original/ --path2 original
	@echo -e "\n\
    **** THE FOLLOWING TEXT SHOULD APPEAR AND IS NOT AN ERROR ****\n\
    ERROR - EXITING, 2: Read of pedigree in file "original/pedfile.dat" failed, in family 9, individual(s) 17 disconnected from the rest of the pedigree"
	- $(KELVIN_ROOT)/kdiff.pl --ped --c1 /dev/null --c2 /dev/null --path1 original/ --path2 pedfile_disconnected
	@echo -e "\n\
    **** THE FOLLOWING TEXT SHOULD APPEAR AND IS NOT AN ERROR ****\n\
    ERROR - EXITING, 2: Read of pedigree in file "original/pedfile.dat" failed, pedfile_extra-eol-stuff/pedfile.dat, line 13: unexpected number of columns, can't guess format"
	- $(KELVIN_ROOT)/kdiff.pl --ped --c1 /dev/null --c2 /dev/null --path1 original/ --path2 pedfile_extra-eol-stuff
	@echo -e "\n\
    **** THE FOLLOWING TEXT SHOULD APPEAR AND IS NOT AN ERROR ****\n\
    2: Ped "1" ind "6" has different values for marker 4 ("Loci_5") - 1:2 1 vs 2:1 1\n\
    2: Ped "2" ind "2" has different values for trait 1 - 1:1 vs 2:2\n\
    2: Ped "3" ind "6" has different values for sex - 1:2 vs 2:1\n\
    2: Ped "4" ind "4" has different values for sex - 1:1 vs 2:2\n\
    2: Ped "4" ind "5" has different values for sex - 1:2 vs 2:1\n\
    2: Ped "4" ind "6" has different values for dadid - 1:4 vs 2:5\n\
    2: Ped "4" ind "6" has different values for momid - 1:5 vs 2:4\n\
    2: Ped "9" ind "13" has different values for marker 2 (\"Loci_3\") - 1:1 1 vs 2:1 3\n\
    ERROR - EXITING, Files are different"
	- $(KELVIN_ROOT)/kdiff.pl --ped --c1 /dev/null --c2 /dev/null --path1 original/ --path2 pedfile_many-non-structural-changes
	@echo -e "\n\
    **** THE FOLLOWING TEXT SHOULD APPEAR AND IS NOT AN ERROR ****\n\
    2: Ped "9" ind "13" has different values for marker 2 ("Loci_3") - 1:1 1 vs 2:1 3\n\
    ERROR - EXITING, Files are different"
	- $(KELVIN_ROOT)/kdiff.pl --ped --c1 /dev/null --c2 /dev/null --path1 original/ --path2 pedfile_missing-some-originals
	@echo -e "\n\
    **** THE FOLLOWING TEXT SHOULD APPEAR AND IS NOT AN ERROR ****\n\
    2: Pedigree "5" not found, skipping!\n\
    2: Ped "9" ind "13" has different values for marker 2 ("Loci_3") - 1:1 1 vs 2:1 3\n\
    1: Pedigree "25" not found, skipping!\n\
    ERROR - EXITING, Files are different"
	- $(KELVIN_ROOT)/kdiff.pl --ped --c1 /dev/null --c2 /dev/null --path1 original/ --path2 pedfile_missing_pedigree
	@echo -e "\n\
    **** THE FOLLOWING TEXT SHOULD APPEAR AND IS NOT AN ERROR ****\n\
    ERROR - EXITING, 2: Read of pedigree in file "original/pedfile.dat" failed, in family 2, individual(s) 19 disconnected from the rest of the pedigree"
	- $(KELVIN_ROOT)/kdiff.pl --ped --c1 /dev/null --c2 /dev/null --path1 original/ --path2 pedfile_out-of-order
	$(KELVIN_ROOT)/kdiff.pl --ped --c1 /dev/null --c2 /dev/null --path1 pedfile_oversize-and-string-ids/ --path2 pedfile_oversize-and-string-ids
	@echo -e "\n\
    **** THE FOLLOWING TEXT SHOULD APPEAR AND IS NOT AN ERROR ****\n\
    ERROR - EXITING, 2: Read of pedigree in file "original/pedfile.dat" failed, pedfile_short-rows/pedfile.dat, line 98: too few columns in pedigree file"
	- $(KELVIN_ROOT)/kdiff.pl --ped --c1 /dev/null --c2 /dev/null --path1 original/ --path2 pedfile_short-rows
	@echo -e "\n\
    **** THE FOLLOWING TEXT SHOULD APPEAR AND IS NOT AN ERROR ****\n\
    2: Ped "7" ind "1" has different values for marker 2 ("Loci_3") - 1:2 1 vs 2:T A\n\
    2: Ped "7" ind "2" has different values for marker 2 ("Loci_3") - 1:1 2 vs 2:A T\n\
    2: Ped "7" ind "3" has different values for marker 2 ("Loci_3") - 1:2 1 vs 2:T A\n\
    2: Ped "7" ind "4" has different values for marker 2 ("Loci_3") - 1:2 1 vs 2:T A\n\
    2: Ped "7" ind "5" has different values for marker 2 ("Loci_3") - 1:2 2 vs 2:T T\n\
    2: Ped "7" ind "6" has different values for marker 2 ("Loci_3") - 1:1 1 vs 2:A A\n\
    2: Ped "7" ind "7" has different values for marker 2 ("Loci_3") - 1:1 1 vs 2:A A\n\
    2: Ped "7" ind "8" has different values for marker 2 ("Loci_3") - 1:1 1 vs 2:A A\n\
    2: Ped "7" ind "9" has different values for marker 2 ("Loci_3") - 1:1 1 vs 2:A A\n\
    2: Ped "7" ind "10" has different values for marker 2 ("Loci_3") - 1:1 1 vs 2:A A\n\
    2: Ped "7" ind "11" has different values for marker 2 ("Loci_3") - 1:1 2 vs 2:A T\n\
    2: Ped "7" ind "12" has different values for marker 2 ("Loci_3") - 1:2 1 vs 2:T A\n\
    2: Ped "7" ind "13" has different values for marker 2 ("Loci_3") - 1:2 2 vs 2:T T\n\
    2: Ped "9" ind "13" has different values for marker 2 ("Loci_3") - 1:1 1 vs 2:1 3\n\
    ERROR - EXITING, Files are different"
	- $(KELVIN_ROOT)/kdiff.pl --ped --c1 /dev/null --c2 /dev/null --path1 original/ --path2 pedfile_text-alleles
	@echo -e "\n\
    **** THE FOLLOWING TEXT SHOULD APPEAR AND IS NOT AN ERROR ****\n\
    1: Map doesn't provide phys, will not be compared!\n\
    2: Order of markers in maps is different!\n\
    2: Marker Loci_1 has a different value for chr - 1:40 vs 2:22\n\
    2: Marker Loci_10 has a different value for chr - 1:40 vs 2:22\n\
    1: Marker "Loci_11" not found in map, skipping!\n\
    2: Marker Loci_3 has a different value for chr - 1:40 vs 2:22\n\
    2: Marker Loci_4 has a different value for chr - 1:40 vs 2:22\n\
    2: Marker Loci_5 has a different value for chr - 1:40 vs 2:22\n\
    2: Marker Loci_6 has a different value for chr - 1:40 vs 2:22\n\
    2: Marker Loci_7 has a different value for avgpos - 1:1.99 vs 2:1.98\n\
    2: Marker Loci_7 has a different value for chr - 1:40 vs 2:22\n\
    2: Marker Loci_8 has a different value for chr - 1:40 vs 2:22\n\
    2: Marker Loci_9 has a different value for chr - 1:40 vs 2:22\n\
    ERROR - EXITING, Files are different"
	- $(KELVIN_ROOT)/kdiff.pl --map --c1 /dev/null --c2 /dev/null --path1 original/ --path2 maps_different
	@echo -e "\n\
    **** THE FOLLOWING TEXT SHOULD APPEAR AND IS NOT AN ERROR ****\n\
    BLAH\n\
    ERROR - EXITING, Files are different"
	- $(KELVIN_ROOT)/kdiff.pl --freq --c1 /dev/null --c2 /dev/null --path1 original/ --path2 markers_different

