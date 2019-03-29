setenv testpath ~/kelvin/test-suite/
setenv binarypath ~/kelvin/32/
setenv polynomialDebugFlag 1
setenv test 10; setenv binary kelvin; setenv testsuite 1p1m
mkdir $binarypath$test-$binary-$testsuite/
cd $binarypath$test-$binary-$testsuite/
cp $testpath$testsuite/kelvin.conf .
cp $testpath$testsuite/mapfile.dat .
cp $testpath$testsuite/markers.dat .
cp $testpath$testsuite/pedpost.dat .
cp $testpath$testsuite/datafile.dat .
$binarypath$binary kelvin.conf >& test.log
