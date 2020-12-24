# Kelvin
# Copyright 2010, The Research Institute At Nationwide Children's Hospital
# Permission granted to distribute and use for non-profit educational purposes
# only.

## Directory into which compiled executables and scripts will be installed.
ifndef BINDIR
  BINDIR=/usr/local/bin/kelvin
endif

## Directory on your PATH into which main Kelvin programs will be optionally
## symlinked. (Comment out to not put Kelvin on your PATH)
## Executables/scripts so symlinked: `Kelvin`, `calc_updated_ppl`, and one
## additional script each for exotic operating modes if those are installed.
PATHDIR=/usr/local/bin

## User and group IDs by which installed executables and scripts will be owned.
OWNER=root
GROUP=root

## The C compiler to be used to build executables. Pick one.
## GCC (GNU C Compiler)
# The following work fine
#CC := gcc
#CC := clang 
#CC := gcc-6
#CC := gcc-7
#CC := gcc-8
CC := gcc-9
# The following fail due to lost '.h' files
#CC := gcc-4.9
#CC := gcc-5
# The following fail due to duplicate symbols (which appears to be a real problem)
#CC := gcc-10
# The following is currently untested
## ICC (Intel C Compiler)
# CC := icc

## GCC optimization level, 0=none, 1=default, 2=some (recommended), 3=all
GCCOPT := 2


include Makefile.main
