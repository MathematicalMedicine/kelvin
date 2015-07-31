# Kelvin
# Copyright 2010, The Research Institute At Nationwide Children's Hospital
# Permission granted to distribute and use for non-profit educational purposes
# only.

## Directory into which compiled executables and scripts will be installed.
## For LKS, this should generally NOT be a directory in $PATH.
ifndef BINDIR
  BINDIR=/usr/local/share/kelvin
endif

## Directory on your system's $PATH into which Kelvin-LKS front end scripts 
## will be symlinked
PATHDIR=/usr/local/bin

## User and group IDs by which installed executables and scripts will be owned.
OWNER=root
GROUP=root

## The C compiler to be used to build executables. Pick one.
## GCC (GNU C Compiler)
CC := gcc
## ICC (Intel C Compiler)
# CC := icc

## GCC optimization level, 0=none, 1=default, 2=some (recommended), 3=all
GCCOPT := 2

## LKS Portable Edition site-wide configuration
## These values can be changed after installation by editing
## $(BINDIR)/site_settings.sh
# MySQL distribution location, for setting up transient MySQL servers
LKSPE_MYSQL_BASE="/opt/mysql-5.6"

# Directories that need to also go into $TOOLPATH. This is for if you had more
# than one directory in DEFAULT_TOOLPATH in old versions of Portable Edition;
# typically it's only used at sites that also need to include $MYSQL_BASE/bin
LKSPE_EXTENDED_TOOLPATH=""

# qsub command line arguments used to specify that a job will itself submit
# jobs. If left blank, then it's assumed all jobs can submit subsequent jobs.
LKSPE_JOB_SUBMITS_JOBS=$(shell if [[ "$$HOSTNAME" == "Levi-Montalcini"* ]]; then echo "-q headnode -P headnode"; fi)



include Makefile.main
