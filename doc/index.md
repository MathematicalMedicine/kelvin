Introduction to Kelvin
======================

Kelvin is a program suite for analysis of genetic data. It is based on the PPL framework<sup>[1, 2, 3](#references)</sup>, and produces output on the posterior probability (0,..,1) scale.

Kelvin consists of a binary core and Perl front-end combination genetic data analysis program based on the PPL framework, with likelihoods calculated via the Elston-Stewart algorithm<sup>[4](#references)</sup>. It is stable software, and has been tested on a number of hardware platforms.

A discussion of the guiding philosophy of Kelvin and details on the underlying statistical methods can be found in the following reference:
> Vieland VJ, Huang Y, Seok S-C, Burian J, Catalyurek U, O’Connell J, Segre A, Valentine-Cooper W. [KELVIN: a software package for rigorous measurement of statistical evidence in human genetics.](https://doi.org/10.1159/000330634) Hum Hered 72(4):276-288, 2011. PMID: [22189470](https://pubmed.ncbi.nlm.nih.gov/22189470/), [PMC3267994](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3267994/), doi:[10.1159/000330634](https://doi.org/10.1159/000330634)

Details on the usage of Kelvin can be found in our [detailed usage documentation](kelvin-usage-details.html); this document provides a general overview and "getting started" information.

Prerequisites
=============

Kelvin has been tested and run on several platforms, but the reference and development platform is CentOS 6 (or any other Linux distribution of similar vintage).

To install Kelvin, you will also need a working C compiler (GCC will do and is tested; ICC (the Intel C Compiler) has also been tested). You will also almost certainly want libgsl (the GNU Scientific Library); compiling without GSL is an option but not supported by default. pkg-config is also needed; normally this is included in any install of Linux development tools, but we've seen instances where it wasn't.

Running Kelvin requires libgsl (if compiled with same) and Perl 5.8 (or any later version).


Installation
============

1. Edit the Makefile as follows:
    `BINDIR`: This should point to where Kelvin and related modules and utility scripts should be located. The default is `/usr/local/share/kelvin`
    `PATHDIR`: This should point to a directory on your $PATH where the Kelvin program will be linked. The default is `/usr/local/bin`
    `OWNER` and `GROUP`: These should be the owner and group IDs for the Kelvin programs and utility files. The defaults are `root` for both.

2. Run `make install`. Kelvin will be built, assembled, and installed in the location you specified in the Makefile.

3. (optional) Verify the build worked by running `make check` (for a quick check) or `make test` (for a more involved one).

Uninstallation may be done by running `make uninstall`; this simply deletes all files that were installed.


Using Kelvin
============

Kelvin requires four input data files, inspired by the de-facto standard formats employed by the LINKAGE program<sup>[5](#references)</sup>. Details of the formats can be found in our detailed usage documentation under ["Input Data File Formats"](kelvin-usage-details.html#input-data-file-formats). Examples are given here showing an affected sib-pair family with three markers:

* Pedigree File - This contains phenotypic and genotypic information. This will nearly always be in pre-MAKEPED format. (There are some cases where Kelvin will require post-MAKEPED format).

        fam1 papa 0 0 1 1 2 2 1 2 1 1
        fam1 mama 0 0 2 1 1 1 1 2 1 1
        fam1 kid1 papa mama 2 2 1 2 2 2 1 1
        fam1 kid2 papa mama 1 2 2 1 1 1 1 1

* Locus File (also called Data File) - Describes marker column order in the pedigree file, starting with the position of the trait locus.

        T Trait
        M MRK_1
        M MRK_2
        M MRK_3

* Frequency File - Gives the allele frequencies for the markers.

        M MRK_1
        F 0.3 0.7
        M MRK_2
        F 0.35 0.65
        M MRK_3
        F 0.7 0.3

* Map File - Gives the chromosomal position of the markers.

        CHROMOSOME      MARKER          POSITION
        3               MRK_1           0.33
        3               MRK_2           0.66
        3               MRK_2           0.99

Kelvin also requires a _configuration file_. Additional details on creation of same can be found in the detailed usage documentation under ["Preparation and Analysis Considerations"](kelvin-usage-details.html#preparation-and-analysis-considerations) and the ["Configuration File Reference"](kelvin-usage-details.html#configuration-file-reference).

Once your configuration file is created, invoke Kelvin as `Kelvin <configfile>`.


Visualizing Results
===================

Kelvin-formatted PPL output can be easily visualized using our graphing application, Kelviz. Kelviz is distributed separately<!--; information on same can be found in the [Kelviz documentation](http://kelvin.mathmed.org/static/doc/kelviz/index.html) and downloads can be found on the [Kelvin website](https://kelvin.mathmed.org)-->.


References
==========
1. Smith C.A.B. [Some comments on the statistical methods used in linkage investigations.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1932043/) Am J Hum Genet 11(4):289-304, 1959. PMID: [13831904](https://pubmed.ncbi.nlm.nih.gov/13831904/), [PMC1932043](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1932043/)
2. Vieland VJ. [Bayesian linkage analysis, or: How I learned to stop worrying and love the posterior probability of linkage.](https://doi.org/10.1086/302076) Am J Hum Genet 63(4):947-954, 1998. PMID: [9758634](https://pubmed.ncbi.nlm.nih.gov/9758634/), [PMC1377508](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1377508/), doi:[10.1086/302076](https://doi.org/10.1086/302076)
3. Vieland VJ. [Thermometers: Something for statistical geneticists to think about.](https://doi.org/10.1159/000093775) Hum Hered 61:144-156, 2006. PMID: [16770079](https://pubmed.ncbi.nlm.nih.gov/16770079/), doi:[10.1159/000093775](https://doi.org/10.1159/000093775)
4. Elston R.C. & Stewart J. [A general model for the genetic analysis of pedigree data.](https://doi.org/10.1159/000152448) Hum Hered 21:523-542, 1971. PMID: [5149961](https://pubmed.ncbi.nlm.nih.gov/5149961/), doi:[10.1159/000152448](https://doi.org/10.1159/000152448)
5. Ott J. [A computer program for linkage analysis of general human pedigrees.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1685114/) An J Hum Genet 28(5):528-529, 1976. PMID: [984049](https://pubmed.ncbi.nlm.nih.gov/984049/), [PMC1685114](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1685114/)
