#!/bin/csh -f
# compile.sh
set verbose

# Compiler can be passed as the first argument (defaults to gfortran):
#   ./compile.sh           # uses gfortran
#   ./compile.sh pgf90     # uses PGI
 set compiler = gfortran
 if ($#argv >= 1) then
    set compiler = "$argv[1]"
 endif
 if ($compiler == pgf90) then
    set flags = "-g -C -byteswapio -Ktrap=fp -O0"
 else if ($compiler == gfortran) then
    set flags = "-g -fallow-argument-mismatch"
 else
    echo "*** bad compiler = " $compiler
    exit
 endif

#/bin/rm mosaic.x

cd  compile
/bin/rm  *.f90  Makefile*

cp -p  ../datamodules/*.f90  .
cp -p  ../main/*.f90  .
cp -p  ../gas/*.f90  .
cp -p  ../solver/*.f  .
cp -p  ../aerosol/*.f90  .
cp -p  ../cloud/*.f90  . 

#BSINGH mosaic_support.f90 has to go through C pre-processor (cpp)
/bin/rm module_mosaic_support.f90
#cpp ../main/module_mosaic_support.f90 module_mosaic_support.f90
cpp -P -traditional-cpp ../main/module_mosaic_support.f90 module_mosaic_support.f90

cp -p  ../Makefile  .

make mosaic.x  COMP=$compiler  FLAGS="$flags" --always-make


make libmosaic.a COMP=$compiler  FLAGS="$flags" --always-make
#/bin/mv  mosaic.x  ..

# Export build for linking against MOSAIC (e.g. PartMC)
mkdir -p ../include ../lib
cp -p *.mod          ../include/
cp -p libmosaic.a    ../lib/

cd  ..

unset verbose
