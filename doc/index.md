Introduction to Kelvin
======================

Kelvin is a program suite for analysis of genetic data. It is based on the PPL framework<sup>[1, 2, 3](#references)</sup>, and produces output on the posterior probability (0,..,1) scale.

Kelvin comes in two primary forms:

* The original Kelvin program (hereafter "Kelvin-original") - a binary core and Perl front-end combination genetic data analysis program based on the PPL framework, with likelihoods calculated via the Elston-Stewart algorithm<sup>[4](#references)</sup>. The bulk of our available documentation covers this implementation.
* Kelvin and Likelihood Server, Portable Edition (hereafter "Kelvin-LKS") - the original Kelvin program, paired with additional analysis programs by way of a "likelihood server" (implemented using MySQL) to enable alternative likelihood calculation algorithms (MC-MC using blocked Gibbs sampling<sup>[5](#references)</sup> derived from [JPSGCS](http://balance.med.utah.edu/wiki/index.php/JPSGCS) and Lander-Green<sup>[6](#references)</sup> via Merlin<sup>[7](#references)</sup>) and parallelization of analysis.

Kelvin-original is stable software, and has been tested on a number of hardware platforms. Kelvin-LKS is experimental alpha-quality software and currently only supports linkage analysis (linkage disequilibrium analysis is presently being implemented). It's under active development, and we endeavor to provide prompt support, but it can and will break under various circumstances.

A discussion of the guiding philosophy of Kelvin and details on the underlying statistical methods can be found in the following reference:
> [Vieland, V.J., et al. Kelvin: a Software Package for Rigorous Measurement of Statistical Evidence in Human Genetics. _Hum Hered_ 2011;72(4):276-88. Epub 2011 Dec 23. PMID:22189470](http://kelvin.mathmed.org/static/Kelvin.pdf)



Prerequisites
=============

Kelvin-original
---------------

Kelvin-original has been tested and run on several platforms, but the reference and development platform is CentOS 6 (or any other Linux distribution of similar vintage).

To install Kelvin, you will also need a working C compiler (GCC will do and is tested; ICC (the Intel C Compiler) has also been tested). You will also almost certainly want libgsl (the GNU Scientific Library); compiling without GSL is an option but not supported by default. pkg-config is also needed; normally this is included in any install of Linux development tools, but we've seen instances where it wasn't.

Running Kelvin requires libgsl (if compiled with same) and Perl 5.8 (or any later version)

Kelvin-LKS
----------

Kelvin-LKS, in addition to the requirements for Kelvin-original, requires an Open Grid Scheduler (or Sun Grid Engine, or another descendant thereof) cluster. Furthermore, the following is necessary of the cluster:

* Some jobs need to be able to submit subsequent jobs. Kelvin-LKS by default will assume all nodes can do this, but also supports configurations in which additional options to qsub are necessary to enable that job to submit other jobs.
* Some nodes will have to be made available as "database" nodes - these are nodes that will temporarily run MySQL during analysis (MySQL does not need to be preinstalled on the nodes). These have to be allocated and tracked by the scheduler.

Build requirements are:

* GCC and G++
* pkg-config (should be present, but we've seen platforms where it isn't)
* libgsl, with headers (unlike Kelvin-original, it is not optional for Kelvin-LKS)
* libmysqlclient, with headers
* an Internet connection (so that we can download and build Merlin)

Additional software requirements include:

* Java 6 (or a later version)
* MySQL client version 8 and associated libraries
* Perl 5.8 (or a later version), with the following modules:
    * DBI
    * DBD::mysql
    * DBD::SQLite
    * List::MoreUtils

Database nodes will need a prebuilt binary distribution of MySQL Server version 8, NOT SET UP. This can be downloaded from mysql.com (look for "MySQL Community Server"'s "Linux - Generic", "Compressed TAR Archive").


Installation
============

For Kelvin-LKS and Kelvin-original together, installation requires several steps:

1. Select which nodes in your cluster are to be "database" nodes. Unpack your MySQL prebuilt distribution into the same directory on each node.

2. Add a "database" INT resource to the scheduler and add it to each DB node's complex values. A Perl script (`LKS_setupSGEDB.pm`) is provided that can do this for you; run it with the `--help` option for guidance.

3. Edit the Makefile as follows:
    `BINDIR`: This should point to where Kelvin and related modules and utility scripts should be located. The default is `/usr/local/share/kelvin`
    `PATHDIR`: This should point to a directory on your $PATH where the Kelvin program will be linked. The default is `/usr/local/bin`
    `OWNER` and `GROUP`: These should be the owner and group IDs for the Kelvin programs and utility files. The defaults are `root` for both.
    `LKSPE_MYSQL_BASE`: This is the directory where the MySQL Server binary distribution is located on each "database" node.
    `LKSPE_JOB_SUBMITS_JOBS`: This should be set if your cluster requires some sort of additional command-line option(s) for qsub to indicate jobs that submit other jobs.

4. Run `make install-lks`. Kelvin-LKS, Kelvin-original, and associated programs will be built, assembled, and installed in the location you specified in the Makefile.

5. (optional) Verify that setup worked by running one or both of the two Acid Tests - these are preassembled analyses stored in tarballs in `PATHDIR` named `sa_dt-acid-test.tar.gz` and `merlin-only-sadt-acid-test.tar.gz`. The former tests the whole thing; the latter verifies that Merlin integration is working.

Uninstallation may be done by running `make uninstall-lks`; this simply deletes all files that were installed.


If you are only interested in Kelvin-original, installation is simpler:

1. Edit the Makefile as per Step 3 above, ignoring those variables that start with `LKSPE_`.

2. Run `make install`. Kelvin will be built, assembled, and installed in the location you specified in the Makefile.

3. (optional) Verify the build worked by running `make check` (for a quick check) or `make test` (for a more involved one).

Uninstallation may be done by running `make uninstall`; this simply deletes all files that were installed.


Input Files for Kelvin
======================

Kelvin - in both forms - requires four input data files, inspired by the de-facto standard formats employed by the LINKAGE program<sup>[8](#references)</sup>. Examples are given here showing an affected sib-pair family with three markers:

* Pedigree File - This contains phenotypic and genotypic information. This will nearly always be in pre-MAKEPED format. (There are some cases where Kelvin-original will require post-MAKEPED format).

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

Kelvin-original also requires a _configuration file_. Additional details on creation of same and the use of Kelvin-original can be found in [the Kelvin-original usage documentation](kelvin-original-usage.html).

Kelvin-LKS may also sometimes require a _loop-breaker file_. This consists of a list of individuals (specified via family id followed by individual id) that indicate where loops in pedigrees should be broken. This is necessary for any pedigrees with loops. **Loops should NOT be broken in advance!**


Using Kelvin
============

Kelvin-original, in addition to its input files, requires a _configuration file_ for usage. Once the configuration file is created, it is invoked as `Kelvin <configfile>`. Additional details on creation of same and the use of Kelvin-original can be found in [the Kelvin-original usage documentation](kelvin-original-usage.html).

Kelvin-LKS has a multi-step startup and invocation process:

* Create a working directory to run your analysis in that is accessible from all nodes on your cluster.
* (optional) Copy or move your input data files to this working directory.
* In the working directory, run `ready_kelvin_lks_analysis.sh`. This will create two new scripts in the folder - `control.sh` and `settings.sh` - and create a `processing_output` directory.
* Edit parameters in `settings.sh`. Some are necessary (`CHRS` and the `INPUT_*` parameters are there to indicate where and how to find your data), some are helpful (`OPERATOR`, `DEVELOPER`, and `ANALYSIS` will allow you to receive helpful progress emails and warn you when things go wrong), and some are strictly optional (such as the `MCMC_*` parameters. They're all documented right there in the file, so that should help.
* Run `control.sh`.

Final results of a Kelvin-LKS analysis will be in the file `pooled/pooled.ppl.out`.


Visualizing Results
===================

Kelvin-formatted PPL output can be easily visualized using our graphing application, Kelviz. Kelviz is distributed separately; information on same can be found [in the Kelviz documentation](http://kelvin.mathmed.org/static/doc/kelviz/index.html) and downloads can be found [on the Kelvin website](http://kelvin.mathmed.org).


References
==========
1. Smith, C.A.B. Some comments on the statistical methods used in linkage investigations. _Am J Hum Genet_ 11, 289-304 (1959).
2. Vieland, V.J. Bayesian linkage analysis, or: how I learned to stop worrying and love the posterior probability of linkage. _Am J Hum Genet_ 63, 947-54 PMID: 9758634 (1998).
3. Vieland, V.J. Thermometers: something for statistical geneticists to think about. _Hum Hered_ 61, 144-56 PMID: 16770079 (2006).
4. Elston, R.C. & Stewart, J. A general model for the genetic analysis of pedigree data. _Hum Hered_ 21, 523-42 (1971).
5. Thomas A., Gutin A., Abkevich V., and Bansal A. (2000). Multilocus linkage analysis by blocked Gibbs sampling. _Stat. Comput._ 10, 259-269.
6. Lander, E.S. & Green, P. Construction of multilocus genetic linkage maps in humans. _Proc Natl Acad Sci U S A_ 84, 2363-7 (1987). 
7. Abecasis GR, Cherny SS, Cookson WO and Cardon LR (2002). Merlin-rapid analysis of dense genetic maps using sparse gene flow trees. _Nat Genet_ 30:97-101.
8. Ott, J. (1976). A computer program for linkage analysis of general human pedigrees. _Am. J. Hum. Genet._ 28, 528-529.
