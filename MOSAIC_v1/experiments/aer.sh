# compile.sh
set verbose

/msrc/lapps/pgi-7.1/linux86-64/7.1/bin/pgf90 -o aer_input.x -g -C -Mextend=132 aer_input.f90

aer_input.x

unset verbose
