# this program is written using MacOS as operative system. Here LAPACK is integrated in the OS 
# and to enable it it is sufficient to add the option -framework Accelerate

gfortran Ex05-MatteoBortoletto.f90 -o Ex05.x -O3 -framework Accelerate

./Ex05.x