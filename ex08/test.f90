program test_density_matr

    use Utilities
    use debugger 
    use DensityMatrices

    implicit none 

    integer :: N, D, info
    complex*16, dimension(:), allocatable :: psi
    complex*16, dimension(:,:), allocatable :: rho, rho_A, rho_B
    complex*16 :: trace

    ! ---- set N=2 and D=2 ------------------------------------------------------
    N = 2
    D = 2
    ! ---------------------------------------------------------------------------

    ! ---- allocate memory ------------------------------------------------------
    allocate(psi(D**N), stat=info)
    call checkpoint(debug=info/=0, message="Error in memory allocation.", &
                    end_program=.true.)
    ! ---------------------------------------------------------------------------

    ! ---- initialize the state in case of separable state ----------------------
    psi(1) = 1.d0/sqrt(2.d0) 
    psi(2) = 1.d0/sqrt(2.d0)            
    psi(3) = 0.d0
    psi(4) = 0.d0
    ! ---------------------------------------------------------------------------

    ! ---- compute the density matrices -----------------------------------------
    allocate(rho(D**N, D**N), stat=info)
    call checkpoint(debug=info/=0, message="Error in memory allocation.", &
                    end_program=.true.)
    call OuterProd(psi, conjg(psi), rho)
    ! ---------------------------------------------------------------------------

    ! ---- check if the trace is 1 ----------------------------------------------
    trace = ComplexSquareMatrixTrace(rho)
    print *, "Trace(rho_sep) =", trace
    ! ---------------------------------------------------------------------------

    ! ---- compute the reduced density matrix for the general density matrix ----
    allocate(rho_A(D**(N-1), D**(N-1)), rho_B(D**(N-1), D**(N-1)), stat=info)
    call checkpoint(debug=info/=0, message="Error in memory allocation.", &
                    end_program=.true.)
    call ReducedDensityMatrixA(rho, D, rho_A)
    call ReducedDensityMatrixB(rho, D, rho_B) 
    ! ---------------------------------------------------------------------------

    ! ---- check if the trace is 1 ----------------------------------------------
    call checkpoint(debug=abs(ComplexSquareMatrixTrace(rho_A)-1).ge.1e-04, &
                    message="Trace(rho_get) is not 1!", &
                    end_program=.false.)
    call checkpoint(debug=abs(ComplexSquareMatrixTrace(rho_B)-1).ge.1e-04, &
                    message="Trace(rho_get) is not 1!", &
                    end_program=.false.)
    print *, "Trace(rho_A) =", ComplexSquareMatrixTrace(rho_A)
    print *, "Trace(rho_B) =", ComplexSquareMatrixTrace(rho_B)
    ! ---------------------------------------------------------------------------

    ! ---- save results ---------------------------------------------------------
    call WriteComplex16Matr(rho, trace, 'rho.txt') 
    call WriteComplex16Matr(rho_A, ComplexSquareMatrixTrace(rho_A), 'rhoA.txt') 
    call WriteComplex16Matr(rho_B, ComplexSquareMatrixTrace(rho_B), 'rhoB.txt') 
    ! ---------------------------------------------------------------------------

    ! ---- deallocate memory ----------------------------------------------------
    deallocate(psi, rho, rho_A, rho_B)
    ! ---------------------------------------------------------------------------

end program test_density_matr