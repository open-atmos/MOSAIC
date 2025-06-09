# compile.sh
set verbose

#/msrc/lapps/pgi-7.1/linux86-64/7.1/bin/pgf90 -o smps.x -g -C -Mextend=132 -r8 smps_input.f90

pgf90 -o smps.x -g -C -Mextend=132 -r8 smps_input.f90

./smps.x

unset verbose
