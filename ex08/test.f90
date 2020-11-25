program test_density_matr

    use Utilities
    use debugger 

    implicit none 


    integer :: N, D
    complex*16, dimension(:), allocatable :: psi
    complex*16, dimension(:,:), allocatable :: rho, rho_A, rho_B
    complex*16 :: trace


    N = 2
    D = 2

    ! ---- allocate memory ----
    allocate(psi(D**N))
    ! -------------------------

    ! ---- initialize the state in case of separable state ----
    psi(1) = 1.d0/sqrt(2.d0) 
    psi(2) = 1.d0/sqrt(2.d0)            
    psi(3) = 0.d0
    psi(4) = 0.d0
    ! ---------------------

    ! ---- compute the density matrices ----
    allocate(rho(D**N, D**N))
    call OuterProd(psi, conjg(psi), rho)
    ! ----

    ! ---- check if the trace is 1 -------
    trace = ComplexSquareMatrixTrace(rho)
    print *, "Trace(rho_sep) =", trace
    ! ------------------------------------

    ! ---- compute the reduced density matrix for the general density matrix ----
    ! Notation:  
    ! rho_gen = rho_{AB} with A = left system, B = right system 
    allocate(rho_A(D**(N-1), D**(N-1)))
    allocate(rho_B(D**(N-1), D**(N-1)))
    call ReducedDensityMatrixA(rho, D, rho_A)
    call ReducedDensityMatrixB(rho, D, rho_B) 

    print *, "Trace(rho_A) =", ComplexSquareMatrixTrace(rho_A)
    print *, "Trace(rho_B) =", ComplexSquareMatrixTrace(rho_B)

    call WriteComplex16Matr(rho, trace, 'rho.txt') 
    call WriteComplex16Matr(rho_A, ComplexSquareMatrixTrace(rho_A), 'rhoA.txt') 
    

end program test_density_matr