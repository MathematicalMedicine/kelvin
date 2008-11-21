#!/bin/bash
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
# Takes a single parameter, which is the .h file of a polynomial to be
# compiled and linked into DLs.
#
name=${1##*/}
name=${name%\.*}

# Level 0 optimization compiles several times faster than level 1, but
# level 1-compiled DLs execute twice as fast as level 0. Level 2 is an
# incredibly slow compile that doesn't perform better than level 1.
optFlag=1

echo Processing ${name} with optimization level ${optFlag}
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
	gcc -g -c -fPIC -shared -I/home/whv001/kelvin/trunk/include -O${optFlag} -o ${src}.o ${src}.c >& ${src}.out
	if test ! -e ${src}.o ; then
	    echo Compile failed for some unknown reason, notifying Bill
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
gcc -g -O${optFlag} -fPIC -shared -o ${name}.so ${name}.o >& ${name}-link.out
if test ! -x ${name}.so ; then
    echo Link of root DL failed for some unknown reason, notifying Bill
    mail -s "Link for root DL ${name} failed on ${HOSTNAME} for some unknown reason" whv001@ccri.net < /dev/null
    exit
fi
mv ${name}-link.out compiled
rm ${name}.o
for dl in ${name}_[0-9]*00.o ; do [ -f $dl ] || continue ;
    dl=${dl##*/}
    dl=${dl%00\.o}
    gcc -g -O${optFlag} -fPIC -shared -o ${dl}00.so ${dl}[0-9][0-9].o >& ${dl}00-link.out
    if test ! -x ${name}.so ; then
	echo Link of branch DL failed for some unknown reason, notifying Bill
	mail -s "Link for branch DL ${dl} failed on ${HOSTNAME} for some unknown reason" whv001@ccri.net < /dev/null
	exit
    fi
    mv ${dl}00-link.out compiled
    rm ${dl}[0-9][0-9].o
done
mv ${name}.h compiled
rm ${name}.linking

