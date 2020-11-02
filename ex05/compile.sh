# the option -framework Accelerate tells to use LAPACK on MacOS
gfortran -c Ex05-MatteoBortoletto.f90 -framework Accelerate

./Ex05.x 