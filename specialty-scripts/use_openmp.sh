#!/bin/bash

# Get all of the GNU C compiler #defines as environment variables
temp=`mktemp 2>/dev/null || mktemp -t 'use_openmp_check'`
echo set +e >$temp
cpp -dM /dev/null | perl -pe "s|#define |export\t|; s|\(||g; s|\)||g; s| |=\"|; s|$|\"|; s|\t| |;" >>$temp
echo set -e >>$temp
source $temp 2>/dev/null
rm $temp

if test "$__VERSION__" \< "4.1"; then
    echo "no"
else
    echo "yes"
fi
