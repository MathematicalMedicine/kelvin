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
make $* CFLAGS=" $ADD_CFLAGS -Wall $WERROR -DGCCOPT=2 -O2 -g -D_REENTRANT -DMEMGRAPH" ADD_LDFLAGS="-g -ldl $PTMALLOC3" kelvin
mv kelvin-$VERSION kelvin-$VERSION-no_GSL

# Set OMP_NUM_THREADS=<something big> for best performance after compiliation of DLs, doesn't do compilation if DL not found.
make clean
make $* CFLAGS=" $ADD_CFLAGS -Wall $WERROR -DGCCOPT=2 -O2 -g -D_REENTRANT $OPENMP -DMEMGRAPH -DUSE_GSL -DPOLYUSE_DL" ADD_LDFLAGS="-g -ldl $PTMALLOC3 -lgsl -lgslcblas -lm" kelvin-$VERSION
mv kelvin-$VERSION kelvin-$VERSION-POLYUSE_DL

# Best reliable method for preparing code for DLs. Need to compile separately and then evaluate.
# Notice no OpenMP, as it gains no speed and loses us lots of memory.
#make clean
#make $* CFLAGS=" $ADD_CFLAGS -Wall $WERROR -DGCCOPT=2 -O2 -g -D_REENTRANT -DMEMGRAPH -DUSE_GSL -DPOLYSTATISTICS -DPOLYUSE_DL -DPOLYCODE_DL -DFAKEEVALUATE  -I/usr/local/mysql/include -I/usr/include/mysql" ADD_LDFLAGS="-g -ldl $PTMALLOC3 -lgsl -lgslcblas -lm  -lklvndb -lmysqlclient -L/usr/local/mysql/lib -L/usr/lib64/mysql"  kelvin
#make $* CFLAGS=" $ADD_CFLAGS -Wall $WERROR -DGCCOPT=2 -O2 -g -D_REENTRANT -DMEMGRAPH -DUSE_GSL -DPOLYSTATISTICS -DPOLYUSE_DL -DPOLYCODE_DL -DFAKEEVALUATE"  ADD_LDFLAGS="-g -ldl $PTMALLOC3 -lgsl -lgslcblas -lm "  kelvin
#mv kelvin-$VERSION kelvin-$VERSION-POLYCODE_DL-FAKEEVALUATE

# Best experimental method for preparing LARGE DLs. Need to compile separately then evaluate.
# Notice no OpenMP, as it gains no speed and loses us lots of memory.
#make clean
#make $* CFLAGS=" $ADD_CFLAGS -Wall $WERROR -DGCCOPT=2 -O2 -g -D_REENTRANT -DMEMGRAPH -DUSE_GSL -DPOLYSTATISTICS -DUSE_SSD -DPOLYUSE_DL -DPOLYCODE_DL -DFAKEEVALUATE"  ADD_LDFLAGS="-g -ldl $PTMALLOC3 -lgsl -lgslcblas -lm "  kelvin
#mv kelvin-$VERSION kelvin-$VERSION-POLYCODE_DL-FAKEEVALUATE-SSD

# Build, code and then compile and evaluate DLs. Works fine for small polynomials.
# Notice no OpenMP, as it gains no speed and loses us lots of memory.
make clean
make $* CFLAGS=" $ADD_CFLAGS -Wall $WERROR -DGCCOPT=2 -O2 -g -D_REENTRANT -DMEMGRAPH -DUSE_GSL -DPOLYSTATISTICS -DPOLYUSE_DL -DPOLYCODE_DL -DPOLYCOMP_DL" ADD_LDFLAGS="-g -ldl $PTMALLOC3 -lgsl -lgslcblas -lm" kelvin
mv kelvin-$VERSION kelvin-$VERSION-POLYCOMP_DL

# Normal every-day use
make clean
make $* CFLAGS=" $ADD_CFLAGS -Wall $WERROR -DGCCOPT=2 -O2 -g -D_REENTRANT $OPENMP -DMEMGRAPH -DUSE_GSL" ADD_LDFLAGS="-g -ldl $PTMALLOC3 -lgsl -lgslcblas -lm" kelvin
cp kelvin-$VERSION kelvin-$VERSION-normal

if test "$HOSTNAME" = "Levi-Montalcini" ; then
    # Likelihood server build.
    make clean
    make $* CFLAGS=" $ADD_CFLAGS -Wall $WERROR -DGCCOPT=2 -O2 -g -D_REENTRANT -DUSE_GSL -DUSE_STUDYDB -DSTUDYDB -I/usr/local/mysql/include -I/usr/include/mysql" ADD_LDFLAGS="-g -ldl $PTMALLOC3 -lgsl -lgslcblas -lm -lrt -lklvndb -lmysqlclient -L/usr/local/mysql/lib -L/usr/lib64/mysql" kelvin
    mv kelvin-$VERSION kelvin-$VERSION-study
fi

make seq_update/calc_updated_ppl
