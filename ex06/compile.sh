gfortran debugger.f90 util.f90 harm_osc_1D.f90 Ex06-MatteoBortoletto.f90 -o Ex06.x -O3 -framework Accelerate
gfortran th_eigenf.f90 -o th_eigenf.x 

./Ex06.x 
./th_eigenf.x

