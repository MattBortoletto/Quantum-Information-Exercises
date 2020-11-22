gfortran -framework Accelerate -L/usr/local/lib -lfftw3 -lm time_evol.f90 debugger.f90 harm_osc_1D.f90 util.f90 ex07t.f90 -o Ex07.x

./Ex07.x