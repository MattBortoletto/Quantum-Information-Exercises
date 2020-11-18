gfortran -L /usr/local/lib -framework Accelerate -lfftw3 -lm -o Ex07.x time_evol.f90 debugger.f90 harm_osc_1D.f90 util.f90 ex07.f90 

./Ex07.x 