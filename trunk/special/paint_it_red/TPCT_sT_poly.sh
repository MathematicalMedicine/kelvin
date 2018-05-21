#!/bin/bash -x
tmpdir=$(mktemp --directory)
olddir=$(pwd)
cd ${tmpdir}
cp ${olddir}/* .
make test
cd $OLDPWD
rm -rf ${tmpdir}
