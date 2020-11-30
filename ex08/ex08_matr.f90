program density_matrices

    use Utilities
    use debugger 
    use DensityMatrices 

    implicit none 

    ! N: number of subsystems
    ! D: dimension of each Hilbert space
    ! info: stat check
    integer :: N, D, info
    ! psi_pure: psi for a separable state
    ! psi_gen: general psi 
    complex*16, dimension(:), allocatable :: psi_sep, psi_gen 
    ! rho_sep: density matrix for a separable state
    ! rho_gen: density matrix for a general state
    ! sep_coeff: coefficients for the separable state
    complex*16, dimension(:,:), allocatable :: rho_sep, rho_gen, sep_coeff, rho_A, rho_B
    ! trace_sep: trace of the density matrix for separable states
    ! trace_gen: trace of the density matrix for general states
    complex*16 :: trace_sep, trace_gen 

    ! ---- ask the user to enter D and N ----------------------------------------------------
    print *, &
    "Please enter the dimension of each Hilbert space and the number of subsystems [D,N]:"
    read *, D, N 
    ! N = 2
    ! D = 3
    ! ---------------------------------------------------------------------------------------

    ! ---- allocate memory ------------------------------------------------------------------
    allocate(psi_sep(D**N), psi_gen(D**N), sep_coeff(D, N), stat=info)
    call checkpoint(debug=info/=0, message="Error in memory allocation.", end_program=.true.)
    ! ---------------------------------------------------------------------------------------

    ! ---- initialize the state in case of separable state ----------------------------------
    call InitPureSepState(sep_coeff) 
    ! --------------------------------------------------------------------------------------

    ! ---- initialize the state in case of not-separable state -----------------------------
    call InitPureGenState(psi_gen)
    ! --------------------------------------------------------------------------------------

    ! ---- initialize the separable state using the coefficients matrix ----
    call BuildSeparableStateFromMatr(sep_coeff, psi_sep)
    ! --------------------------------------------------------------------------------------

    ! ---- compute the density matrices ----------------------------------------------------
    allocate(rho_sep(D**N, D**N))
    allocate(rho_gen(D**N, D**N))
    call OuterProd(psi_sep, conjg(psi_sep), rho_sep)
    call OuterProd(psi_gen, conjg(psi_gen), rho_gen) 
    ! --------------------------------------------------------------------------------------

    ! ---- check if the trace is 1 ---------------------------------------------------------
    trace_sep = ComplexSquareMatrixTrace(rho_sep)
    trace_gen = ComplexSquareMatrixTrace(rho_gen)
    print *, "Trace(rho_sep) =", trace_sep
    print *, "Trace(rho_gen) =", trace_gen
    call checkpoint(debug=abs(trace_sep-1).ge.1e-04, message="Trace(rho_sep) is not 1!", &
                    end_program=.false.)
    call checkpoint(debug=abs(trace_gen-1).ge.1e-04, message="Trace(rho_get) is not 1!", &
                    end_program=.false.)
    ! --------------------------------------------------------------------------------------

    ! ---- compute the reduced density matrix for the general density matrix ---------------
    ! Notation: rho_gen = rho_{AB} with A = left system, B = right system 
    allocate(rho_A(D**(N-1), D**(N-1)))
    allocate(rho_B(D**(N-1), D**(N-1)))
    call ReducedDensityMatrixA(rho_gen, D, rho_A)
    call ReducedDensityMatrixB(rho_gen, D, rho_B) 
    ! --------------------------------------------------------------------------------------

    ! ---- check if the trace is 1 ---------------------------------------------------------
    print *, "Trace(rho_A) =", ComplexSquareMatrixTrace(rho_A)
    print *, "Trace(rho_B) =", ComplexSquareMatrixTrace(rho_B)
    ! --------------------------------------------------------------------------------------

    ! ---- deallocate memory ---------------------------------------------------------------
    deallocate(psi_sep, psi_gen, sep_coeff, rho_sep, rho_gen, rho_A, rho_B) 
    ! --------------------------------------------------------------------------------------

end program density_matrices