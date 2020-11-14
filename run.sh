#!/bin/sh

#---------------------------------------------------------------------
# Fourier transform of a real or complex function stored
#   as two column data in a file using FFTW library.
#
# RUN COMMAND:
#     ./spectra.exe <file_name> [<direction>] [<resol_factor>]
#   where:
#     <file_name> - two-column data file with real argument
#       and real or complex value.
#     <direction> - optional transform direction:
#       'f' - forward, 'b' - backward. Forward is used by default.
#     <resol_factor> - an optional integer factor to change result
#       spectra resolution (default value resol_factor = 1).
#       The time duration of the field loaded is incresed up
#       to Tf = T*<resol_factor>. If <resol_factor> > 1
#       the field is set to zero for time after T. Since
#       frequency step of the result spectra is 2Pi/Tf, the
#       spectra resolution multiplies by factor <resol_factor>.
#       If <resol_factor> < 1 it is set to 1.
#---------------------------------------------------------------------
g++ -fopenmp `pkg-config --cflags fftw3` -o spectra.exe \
  spectra.cpp `pkg-config --libs fftw3`;
./spectra.exe test.dat b 3;
./spectra.exe dft_backward_test.dat f;
./spectra.exe test.dat;
./spectra.exe dft_forward_test.dat b 2;
rm *.exe;
