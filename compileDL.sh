#!/bin/bash
set -v
#
# Copyright 2008, Nationwide Children's Research Institute.  All
# rights reserved.  Permission is hereby given to use this software
# for non-profit educational purposes only.
#
# Compile and link the source to a dynamic library, and then remove it.
# Use some flag files to indicate the current status of compiles and
# linkage, since we'll want to distribute the compilation as widely as
# possible.
#
# Takes two parameters:
#
# 1: the name of a polynomial to be compiled and linked into DLs. The
# name will be stripped of any preceeding path and trailing file type,
# i.e. "./CL7_P9C1M_9_T_10_11.h" will be reduced to CL7_P9C1M_9_T_10_11.
#
# 2: the optimization level to use, typically 0 or 1, defaults to 0 if
# not specified.
#
# Presumes the existance of a KELVIN_ROOT environment variable which
# contains the path to the kelvin source hierarchy.
#
# Since the dynamic libraries are created in the default directory,
# LD_LIBRARY_PATH (and DYLD_LIBRARY_PATH under Darwin) should be set as:
#
# export LD_LIBRARY_PATH=.:$KELVIN_ROOT/lib/
#
# If GSL is required and not installed in a normal location, create
# environment variables INCDIR and LIBDIR as paths to the include
# and library locations.
#
if test -z "${INCDIR}" ; then
    INCDIR=/usr/local/include
    LIBDIR=/usr/local/lib
fi

name=${1##*/}
name=${name%\.*}

# Level 0 optimization compiles several times faster than level 1, but
# level 1-compiled DLs execute twice as fast as level 0. Level 2 is an
# incredibly slow compile that doesn't perform better than level 1.
optLevel=${2}
if test -z "${optLevel}" ; then
    optLevel=0
fi

if test "`uname`" = "Darwin" ; then
    shareOpt="-dynamiclib"
else
    shareOpt="-shared"
fi

echo Processing ${name} with optimization level ${optLevel}
if test ! -e compiled ; then
    mkdir compiled
fi

others=""
for src in ${name}{\.,\_[0-9]*}c ; do [ -f $src ] || continue ;
    echo Looking at ${src}
    src=${src##*/}
    src=${src%\.*}
    if test ! -e ${src}.compiling ; then
	touch ${src}.compiling
	echo Compiling ${src}
	gcc -g -c -fPIC -O${optLevel} -I${KELVIN_ROOT}/utils/ -I${INCDIR} -o ${src}.o ${src}.c >& ${src}.out
	if test ! -e ${src}.o ; then
	    echo Compile failed for some unknown reason
	    mail -s "Compile for ${src} failed on ${HOSTNAME} for some unknown reason" whv001@ccri.net < /dev/null
	    exit
	fi
	mv ${src}.c ${src}.out compiled
	rm ${src}.compiling
    else
	others="sure are!"
    fi
done

if test -n "${others}" ; then
    echo Further compiles are in progress, they will have to link
    exit
fi

if test -e ${name}.linking ; then
    echo Already being linked by some other job!
    exit
fi
touch ${name}.linking
gcc -g ${shareOpt} -O${optLevel} -L${KELVIN_ROOT}/lib/ -L${LIBDIR} -o ${name}.so ${name}.o -lgsl -lgslcblas >& ${name}-link.out
if test ! -x ${name}.so ; then
    echo Link of root DL failed for some unknown reason
    mail -s "Link for root DL ${name} failed on ${HOSTNAME} for some unknown reason" whv001@ccri.net < /dev/null
    exit
fi
mv ${name}-link.out compiled
rm ${name}.o
for dl in ${name}_[0-9]*00.o ; do [ -f $dl ] || continue ;
    dl=${dl##*/}
    dl=${dl%00\.o}
    gcc -g -O${optLevel} -L${KELVIN_ROOT}/lib/ ${shareOpt} -L${LIBDIR} -o ${dl}00.so ${dl}[0-9][0-9].o -lgsl -lgslcblas >& ${dl}00-link.out
    if test ! -x ${name}.so ; then
	echo Link of branch DL failed for some unknown reason
	mail -s "Link for branch DL ${dl} failed on ${HOSTNAME} for some unknown reason" whv001@ccri.net < /dev/null
	exit
    fi
    mv ${dl}00-link.out compiled
    rm ${dl}[0-9][0-9].o
done
mv ${name}.h compiled
rm ${name}.linking

