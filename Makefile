# Kelvin
# Copyright (C) 2006, 2007, 2008, 2009, 2010, 2022 Mathematical Medicine LLC
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version.
# You should have received a copy of the GNU General Public License along
# with this program. If not, see <https://www.gnu.org/licenses/>.

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
CC := gcc
#CC := clang 
#CC := gcc-6
#CC := gcc-7
#CC := gcc-8
#CC := gcc-9
# The following fail due to lost '.h' files in the "#include_next fiasco with Xcode"
#CC := gcc-4.9
#CC := gcc-5
# The following only works after adding -fcommon to the CFLAGS (default changed to -fnocommon in GCC 10)
#CC := gcc-10
# The following is currently untested
## ICC (Intel C Compiler)
# CC := icc

## GCC optimization level, 0=none, 1=default, 2=some (recommended), 3=all
GCCOPT := 2


include Makefile.main
