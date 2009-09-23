make clean
make CFLAGS="$* -Wall -Werror -DGCCOPT=3 -O3 -fopenmp -DMEMGRAPH" ADD_LDFLAGS="-ldl" kelvin
mv kelvin kelvin-normal

# Set OMP_NUM_THREADS=<something big> for best performance after compiliation of DLs.
make clean
make CFLAGS="$* -Wall -Werror -DGCCOPT=3 -O3 -fopenmp -DMEMGRAPH -DPOLYUSE_DL" ADD_LDFLAGS="-ldl" kelvin
mv kelvin kelvin-POLYUSE_DL

# Best reliable method for preparing code for DLs. Need to compile separately and then evaluate.
# Notice no OpenMP, as it gains no speed and loses us lots of memory.
make clean
make CFLAGS="$* -Wall -Werror -DGCCOPT=3 -O3 -DMEMGRAPH -DPOLYSTATISTICS -DPOLYUSE_DL -DPOLYCODE_DL -DFAKEEVALUATE" ADD_LDFLAGS="-ldl" kelvin
mv kelvin kelvin-POLYCODE_DL-FAKEEVALUATE

# Best experimental method for preparing LARGE DLs. Need to compile separately then evaluate.
# Notice no OpenMP, as it gains no speed and loses us lots of memory.
make clean
make CFLAGS="$* -Wall -Werror -DGCCOPT=3 -O3 -DMEMGRAPH -DPOLYSTATISTICS -DUSE_SSD -DPOLYUSE_DL -DPOLYCODE_DL -DFAKEEVALUATE" ADD_LDFLAGS="-ldl" kelvin
mv kelvin kelvin-POLYCODE_DL-FAKEEVALUATE-SSD

# Diagnostic for DLs.
# Notice no OpenMP, as it gains no speed and loses us lots of memory.
make clean
make CFLAGS="$* -Wall -Werror -DGCCOPT=3 -O3 -DMEMGRAPH -DPOLYUSE_DL -DPOLYCODE_DL -DPOLYCOMP_DL -DPOLYCHECK_DL" ADD_LDFLAGS="-ldl" kelvin
mv kelvin kelvin-POLYCHECK_DL

# Build, code and then compile and evaluate DLs. Works fine for small polynomials.
# Notice no OpenMP, as it gains no speed and loses us lots of memory.
make clean
make CFLAGS="$* -Wall -Werror -DGCCOPT=3 -O3 -DMEMGRAPH -DPOLYUSE_DL -DPOLYSTATISTICS -DPOLYCODE_DL -DPOLYCOMP_DL" ADD_LDFLAGS="-ldl" kelvin
mv kelvin kelvin-POLYCOMP_DL

# Build, code and then compile and evaluate DLs. Works fine for small polynomials.
# Notice no OpenMP, as it gains no speed and loses us lots of memory.
make clean
make CFLAGS="$* -Wall -Werror -DGCCOPT=3 -O3 -DMEMGRAPH -DPOLYUSE_DL -DPOLYSTATISTICS -DPOLYCODE_DL -DPOLYCOMP_DL -DFAKEEVALUATE" ADD_LDFLAGS="-ldl" kelvin
mv kelvin kelvin-POLYCOMP_DL-FAKEEVALUATE
