program density_matrices

    use Utilities
    use debugger 

    implicit none 

    ! N: number of subsystems
    ! D: dimension of each Hilbert space
    ! ii: variable to loop 
    ! which_subsystem: 1 if A (left), 2 if B (right)
    integer :: N, D, ii, jj
    ! psi_pure: psi for a separable state
    ! psi_gen: general psi 
    complex*16, dimension(:), allocatable :: psi_sep, psi_gen 
    ! rho_sep: density matrix for a separable state
    ! rho_gen: density matrix for a general state
    ! sep_coeff: coefficients for the separable state
    complex*16, dimension(:,:), allocatable :: rho_sep, rho_gen, sep_coeff, rho_A, rho_B
    ! RePart: real part
    ! ImPart: imaginary part
    real*8 :: RePart, ImPart
    ! trace_sep: trace of the density matrix for separable states
    ! trace_gen: trace of the density matrix for general states
    complex*16 :: trace_sep, trace_gen 

    ! ---- ask the user to enter D and N ----
    ! print *, "Please enter the dimension of each Hilbert space (D) and the number of subsystems (N)"
    ! read *, D, N 
    ! ----

    N = 2
    D = 3

    ! ---- allocate memory ----
    allocate(psi_sep(D**N))
    allocate(psi_gen(D**N))
    allocate(sep_coeff(D, N))
    ! -------------------------

    ! ---- initialize the state in case of separable state ----
    do ii = 1, N
        do jj = 1, D 
            call random_number(RePart)
            call random_number(ImPart)
            sep_coeff(jj, ii) = dcmplx(RePart, ImPart) 
        end do 
        sep_coeff(:, ii) = sep_coeff(:, ii) / sqrt(sum(abs(sep_coeff(:, ii))**2))
        if (abs(sum(abs(sep_coeff(:, ii))**2) - 1) .ge. 1e-4) then 
            print *, "Normalization error!"
            stop 
        end if 
    end do 
    ! ---------------------

    ! ---- initialize the state in case of not-separable state ----
    do ii = 1, size(psi_gen)
        call random_number(RePart)
        call random_number(ImPart)
        psi_gen(ii) = dcmplx(RePart, ImPart) 
    end do 
    psi_gen = psi_gen / sqrt(sum(abs(psi_gen)**2))
    ! ----

    ! ---- initialize the separable state using the coefficients matrix ----
    call BuildSeparableState(sep_coeff, psi_sep)
    ! ----

    ! ---- compute the density matrices ----
    allocate(rho_sep(D**N, D**N))
    allocate(rho_gen(D**N, D**N))
    call OuterProd(psi_sep, conjg(psi_sep), rho_sep)
    call OuterProd(psi_gen, conjg(psi_gen), rho_gen) 
    ! ----

    ! ---- check if the trace is 1 -------
    trace_sep = ComplexSquareMatrixTrace(rho_sep)
    trace_gen = ComplexSquareMatrixTrace(rho_gen)
    print *, "Trace(rho_sep) =", trace_sep
    print *, "Trace(rho_gen) =", trace_gen
    ! METTERE UN CHECKPOINT
    ! ------------------------------------

    ! ---- compute the reduced density matrix for the general density matrix ----
    ! Notation:  
    ! rho_gen = rho_{AB} with A = left system, B = right system 
    allocate(rho_A(D**(N-1), D**(N-1)))
    allocate(rho_B(D**(N-1), D**(N-1)))
    call ReducedDensityMatrixA(rho_gen, D, rho_A)
    call ReducedDensityMatrixB(rho_gen, D, rho_B) 

    print *, "Trace(rho_A) =", ComplexSquareMatrixTrace(rho_A)
    print *, "Trace(rho_B) =", ComplexSquareMatrixTrace(rho_B)
    

end program density_matrices