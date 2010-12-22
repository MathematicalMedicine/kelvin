#!/bin/bash -xe
set -xe # Because it isn't enough for Darwin's cron for us to have it on the shebang

# Get all of the GNU C compiler #defines as environment variables
temp=/tmp/`echo $RANDOM`
echo set +e >$temp
cpp -dM /dev/null | perl -pe "s|#define |export\t|; s|\(||g; s|\)||g; s| |=\"|; s|$|\"|; s|\t| |;" >>$temp
echo set -e >>$temp
source $temp
#rm $temp
printenv | sort

if test "$__VERSION__" \< "4.1"; then
    OPENMP=
else
    OPENMP=-fopenmp
fi

case $HOSTNAME in
    Levi-Montalcini* )
	PTMALLOC3="-lptmalloc3 -lpthread"
#        WERROR=-Werror
        ;;
    * )
	PTMALLOC3=
        WERROR=
	;;
esac

if test -n "$1" ; then
    ADD_CFLAGS=$1
    shift
fi

VERSION="`cat .maj`.`cat .min`.`cat .pat`"

make clean
make USE_GSL=no USE_OPENMP=no $* kelvin
mv kelvin-$VERSION kelvin-$VERSION-no_GSL

# Set OMP_NUM_THREADS=<something big> for best performance after compiliation of DLs, doesn't do compilation if DL not found.
make clean
make $* CFLAGS=" $ADD_CFLAGS -Wall $WERROR -DGCCOPT=2 -O2 -g -D_REENTRANT $OPENMP -DMEMGRAPH -DUSE_GSL -DPOLYUSE_DL" ADD_LDFLAGS="-g -ldl $PTMALLOC3 -lgsl -lgslcblas -lm" kelvin
mv kelvin-$VERSION kelvin-$VERSION-POLYUSE_DL

# Build, code and then compile and evaluate DLs. Works fine for small polynomials.
# Notice no OpenMP, as it gains no speed and loses us lots of memory.
make clean
make $* CFLAGS=" $ADD_CFLAGS -Wall $WERROR -DGCCOPT=2 -O2 -g -D_REENTRANT -DMEMGRAPH -DUSE_GSL -DPOLYSTATISTICS -DPOLYUSE_DL -DPOLYCODE_DL -DPOLYCOMP_DL" ADD_LDFLAGS="-g -ldl $PTMALLOC3 -lgsl -lgslcblas -lm" kelvin
mv kelvin-$VERSION kelvin-$VERSION-POLYCOMP_DL

# Normal every-day use
make clean
make $* kelvin
cp kelvin-$VERSION kelvin-$VERSION-normal

if test "$HOSTNAME" = "Levi-Montalcini" ; then
    # Likelihood server build.
    make clean
    make  USE_STUDYDB=yes $* kelvin
    mv kelvin-$VERSION kelvin-$VERSION-study
fi

make seq_update/calc_updated_ppl
