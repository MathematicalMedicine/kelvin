#!/bin/bash
exec make specialty



# Old script follows...
# (note that none of it will run!)

#!/bin/bash -xe
set -xe # Because it isn't enough for Darwin's cron for us to have it on the shebang

# Get all of the GNU C compiler #defines as environment variables
temp=/tmp/`echo $RANDOM`
echo set +e >$temp
cpp -dM /dev/null | perl -pe "s|#define |export\t|; s|\(||g; s|\)||g; s| |=\"|; s|$|\"|; s|\t| |;" >>$temp
echo set -e >>$temp
source $temp
#rm $temp
set | sort

if test "$__VERSION__" \< "4.1"; then
    USE_OPENMP=no
else
    USE_OPENMP=yes
fi

case $HOSTNAME in
    Levi-Montalcini* )
	USE_PTMALLOC3=yes
#        WERROR=-Werror
        ;;
    * )
	USE_PTMALLOC3=no
        WERROR=
	;;
esac

VERSION="`cat .maj`.`cat .min`.`cat .pat`"

make clean
make USE_GSL=no USE_OPENMP=no USE_PTMALLOC3=$USE_PTMALLOC3 $* ENV_CFLAGS=" $WERROR" ENV_LDFLAGS="" kelvin
mv kelvin-$VERSION kelvin-$VERSION-no_GSL

# Set OMP_NUM_THREADS=<something big> for best performance after compiliation of DLs, doesn't do compilation if DL not found.
make clean
make USE_OPENMP=$USE_OPENMP USE_GSL=yes USE_PTMALLOC3=$USE_PTMALLOC3 $* ENV_CFLAGS=" $WERROR -DMEMGRAPH -DPOLYSTATISTICS -DPOLYUSE_DL" ENV_LDFLAGS=" -ldl" kelvin
mv kelvin-$VERSION kelvin-$VERSION-POLYUSE_DL

# Build, code and then compile and evaluate DLs. Works fine for small polynomials.
# Notice no OpenMP, as it gains no speed and loses us lots of memory.
make clean
make USE_OPENMP=no USE_GSL=yes USE_PTMALLOC3=$USE_PTMALLOC3 $* ENV_CFLAGS=" $WERROR -DMEMGRAPH -DUSE_GSL -DPOLYSTATISTICS -DPOLYUSE_DL -DPOLYCODE_DL -DPOLYCOMP_DL" ENV_LDFLAGS=" -ldl"  kelvin
mv kelvin-$VERSION kelvin-$VERSION-POLYCOMP_DL

# Build and code DLs, pretend to evaluate. Best for making hay out of limited memory.
# Notice no OpenMP, as it gains no speed and loses us lots of memory.
make clean
make USE_OPENMP=no USE_GSL=yes USE_PTMALLOC3=$USE_PTMALLOC3 $* ENV_CFLAGS=" $WERROR -DMEMGRAPH -DUSE_GSL -DPOLYSTATISTICS -DUSE_SSD -DPOLYUSE_DL -DPOLYCODE_DL -DFAKEEVALUATE" ENV_LDFLAGS=" -ldl"  kelvin
mv kelvin-$VERSION kelvin-$VERSION-POLYCODE_DL_FAKEEVALUATE_SSD

# Likelihood server build iff mysql exists and is version 5 or better.
over4=$(mysql --version | grep "Distrib [5-9]" || true)
if test -n "$over4"; then
    make clean USE_STUDYDB=yes
    make USE_STUDYDB=yes USE_PTMALLOC3=$USE_PTMALLOC3 $* ENV_CFLAGS=" $WERROR" ENV_LDFLAGS="" kelvin
    cp kelvin-$VERSION kelvin-study
    mv kelvin-$VERSION kelvin-$VERSION-study
fi

# Normal every-day use
make clean
make all
cp kelvin-$VERSION kelvin-normal
cp kelvin-$VERSION kelvin-$VERSION-normal
