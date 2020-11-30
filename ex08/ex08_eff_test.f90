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
    complex*16, dimension(:,:), allocatable :: sep_coeff !, rho_sep, rho_gen
    ! trace_sep: trace of the density matrix for separable states
    ! trace_gen: trace of the density matrix for general states
    !complex*16 :: trace_sep, trace_gen 
    real*8 :: start_gen, end_gen, start_sep, end_sep, time_sep, time_gen

    D = 2

    open(unit=73, file='cpu_times.txt')
    
    do N = 2, 40

        print *, N, "start"

        ! ---- allocate memory ------------------------------------------------------------------
        allocate(psi_sep(D**N), psi_gen(D**N), sep_coeff(D, N), stat=info)
        call checkpoint(debug=info/=0, message="Error in memory allocation.", end_program=.true.)
        ! ---------------------------------------------------------------------------------------

        ! ---- initialize the state in case of separable state ----------------------------------
        call cpu_time(start_sep)
        call InitPureSepState(sep_coeff) 
        call cpu_time(end_sep)
        time_sep = end_sep - start_sep 
        ! --------------------------------------------------------------------------------------

        ! ---- initialize the state in case of not-separable state -----------------------------
        call cpu_time(start_gen)
        call InitPureGenState(psi_gen)
        call cpu_time(end_gen)
        time_gen = end_gen - start_gen 
        ! --------------------------------------------------------------------------------------

        write(73, *) N, time_sep, time_gen 

        ! ---- initialize the separable state using the coefficients matrix ----
        !call BuildSeparableStateFromMatr(sep_coeff, psi_sep)
        ! --------------------------------------------------------------------------------------

        ! ---- deallocate memory ---------------------------------------------------------------
        deallocate(psi_sep, psi_gen, sep_coeff) 
        ! --------------------------------------------------------------------------------------

        print *, N, "end"

    end do 

end program density_matrices