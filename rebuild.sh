#!/bin/bash -xe
set -xe # Because it isn't enough for Darwin's cron for us to have it on the shebang

case $HOSTNAME in
    testbed*|ygg*|RESG8508* )
        OPENMP=
	PTMALLOC3=
        WERROR=-Werror
	;;
    Deimos|*VDI* )
        OPENMP=
        WERROR=
#	ADD_CFLAGS=-DDISTRIBUTION
	ADD_CFLAGS=
	;;
    * )
        OPENMP=-fopenmp
	PTMALLOC3="-lptmalloc3 -lpthread"
        WERROR=-Werror
	;;
esac

if test -n "$1" ; then
    ADD_CFLAGS=$1
    shift
fi

PLATFORM_NAME=`uname -m`-`uname -s`
PLATFORM_NAME=${PLATFORM_NAME// /}

make clean
make $* CFLAGS=" $ADD_CFLAGS -Wall $WERROR -DGCCOPT=2 -O3 -D_REENTRANT -DMEMGRAPH" ADD_LDFLAGS="-ldl $PTMALLOC3" kelvin
mv bin/kelvin.$PLATFORM_NAME kelvin.$PLATFORM_NAME-no_GSL

# Set OMP_NUM_THREADS=<something big> for best performance after compiliation of DLs.
make clean
make $* CFLAGS=" $ADD_CFLAGS -Wall $WERROR -DGCCOPT=2 -O3 -D_REENTRANT $OPENMP -DMEMGRAPH -DUSE_GSL -DPOLYUSE_DL" ADD_LDFLAGS="-ldl $PTMALLOC3 -lgsl -lgslcblas -lm" kelvin
mv bin/kelvin.$PLATFORM_NAME kelvin.$PLATFORM_NAME-POLYUSE_DL

# Best reliable method for preparing code for DLs. Need to compile separately and then evaluate.
# Notice no OpenMP, as it gains no speed and loses us lots of memory.
make clean
make $* CFLAGS=" $ADD_CFLAGS -Wall $WERROR -DGCCOPT=2 -O3 -D_REENTRANT -DMEMGRAPH -DUSE_GSL -DPOLYSTATISTICS -DPOLYUSE_DL -DPOLYCODE_DL -DFAKEEVALUATE" ADD_LDFLAGS="-ldl $PTMALLOC3 -lgsl -lgslcblas -lm" kelvin
mv bin/kelvin.$PLATFORM_NAME kelvin.$PLATFORM_NAME-POLYCODE_DL-FAKEEVALUATE

# Best experimental method for preparing LARGE DLs. Need to compile separately then evaluate.
# Notice no OpenMP, as it gains no speed and loses us lots of memory.
make clean
make $* CFLAGS=" $ADD_CFLAGS -Wall $WERROR -DGCCOPT=2 -O3 -D_REENTRANT -DMEMGRAPH -DUSE_GSL -DPOLYSTATISTICS -DUSE_SSD -DPOLYUSE_DL -DPOLYCODE_DL -DFAKEEVALUATE" ADD_LDFLAGS="-ldl $PTMALLOC3 -lgsl -lgslcblas -lm" kelvin
mv bin/kelvin.$PLATFORM_NAME kelvin.$PLATFORM_NAME-POLYCODE_DL-FAKEEVALUATE-SSD

# Diagnostic for DLs.
# Notice no OpenMP, as it gains no speed and loses us lots of memory.
make clean
make $* CFLAGS=" $ADD_CFLAGS -Wall $WERROR -DGCCOPT=2 -O3 -D_REENTRANT -DMEMGRAPH -DUSE_GSL -DPOLYUSE_DL -DPOLYCODE_DL -DPOLYCOMP_DL -DPOLYCHECK_DL" ADD_LDFLAGS="-ldl $PTMALLOC3 -lgsl -lgslcblas -lm" kelvin
mv bin/kelvin.$PLATFORM_NAME kelvin.$PLATFORM_NAME-POLYCHECK_DL

# Build, code and then compile and evaluate DLs. Works fine for small polynomials.
# Notice no OpenMP, as it gains no speed and loses us lots of memory.
make clean
make $* CFLAGS=" $ADD_CFLAGS -Wall $WERROR -DGCCOPT=2 -O3 -D_REENTRANT -DMEMGRAPH -DUSE_GSL -DPOLYUSE_DL -DPOLYSTATISTICS -DPOLYCODE_DL -DPOLYCOMP_DL" ADD_LDFLAGS="-ldl $PTMALLOC3 -lgsl -lgslcblas -lm" kelvin
mv bin/kelvin.$PLATFORM_NAME kelvin.$PLATFORM_NAME-POLYCOMP_DL

# Build, code and then compile and evaluate DLs. Works fine for small polynomials.
# Notice no OpenMP, as it gains no speed and loses us lots of memory.
make clean
make $* CFLAGS=" $ADD_CFLAGS -Wall $WERROR -DGCCOPT=2 -O3 -D_REENTRANT -DMEMGRAPH -DUSE_GSL -DPOLYUSE_DL -DPOLYSTATISTICS -DPOLYCODE_DL -DPOLYCOMP_DL -DFAKEEVALUATE" ADD_LDFLAGS="-ldl $PTMALLOC3 -lgsl -lgslcblas -lm" kelvin
mv bin/kelvin.$PLATFORM_NAME kelvin.$PLATFORM_NAME-POLYCOMP_DL-FAKEEVALUATE

make clean
make $* CFLAGS=" $ADD_CFLAGS -Wall $WERROR -DGCCOPT=2 -O3 -D_REENTRANT $OPENMP -DMEMGRAPH -DUSE_GSL" ADD_LDFLAGS="-ldl $PTMALLOC3 -lgsl -lgslcblas -lm" kelvin
cp bin/kelvin.$PLATFORM_NAME kelvin.$PLATFORM_NAME-normal

make seq_update/calc_updated_ppl
