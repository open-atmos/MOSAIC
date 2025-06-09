# compile.sh
set verbose

rm mosaic.x

cd  compile
/bin/rm  *.f90  Makefile*

cp -p  ../datamodules/*.f90  .
cp -p  ../main/*.f90  .
cp -p  ../gas/*.f90  .
cp -p  ../solver/*.f90  .
cp -p  ../aerosol/*.f90  .
cp -p  ../cloud/*.f90  . 

cp -p  ../Makefile  .

make  mosaic.x

/bin/mv  mosaic.x  ..

cd  ..

./mosaic.x

unset verbose
