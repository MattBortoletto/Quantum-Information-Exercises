program time_evolution 

    use HarmonicOscillator1D
    use debugger 
    use Utilities 
    use TimeEvolution 

    implicit none 

    complex*16, dimension(:,:), allocatable :: psi_t
    integer :: L, t_len, N, info, jj 
    real*8 :: dx, omega, m, hbar, T, dt, pi !, t_max
    ! H: hamiltonian
    ! laplacian: discretized laplacian 
    ! harmonic_potential: harmonic potential
    complex*16, dimension(:,:), allocatable :: H, laplacian, harmonic_potential
    ! eig: eigenvalues array
    ! x_grid 
    ! p_grid
    real*8, dimension(:), allocatable :: eig, x_grid, p_grid, mean_t
    real*8, dimension(:,:), allocatable :: prob_t, V_t 

    ! ---- initialize parameters -----------------------------------------------------------------
    L = 500
    dx = 0.01
    N = 2*L + 1
    !t_len = 1000
    !t_max = 1.0
    T = 10
    dt = 0.01
    !dt = t_max / t_len 
    t_len = int( T / dt )
    print *, t_len 
    m = 1
    hbar = 1
    omega = 1
    pi = 4.d0 * datan(1.d0)
    ! --------------------------------------------------------------------------------------------

    ! ---- allocate memory -----------------------------------------------------------------------
    allocate(x_grid(N))
    allocate(p_grid(N))
    allocate(laplacian(N, N)) 
    allocate(harmonic_potential(N, N))
    allocate(eig(N))
    allocate(psi_t(N, t_len+1))
    allocate(prob_t(N, t_len+1))
    allocate(V_t(N, t_len+1))
    allocate(mean_t(t_len))
    ! --------------------------------------------------------------------------------------------

    ! ---- compute the grid points ---------------------------------------------------------------
    ! x grid  
    do jj = 1, N 
        ! x_min = -L*dx
        x_grid(jj) = -L*dx + dx*(jj-1)
    end do
    ! p grid 
    do jj = 1, N/2
        ! dp = (2*pi) / (dx*N)
        p_grid(jj) = (jj-1) * (2*pi)/(dx*N)  
    end do 
    do jj = N/2+1, N 
        p_grid(jj) = (jj-1-N) * (2*pi)/(dx*N) 
    end do 
    ! --------------------------------------------------------------------------------------------

    ! ---- compute the hamiltonian ---------------------------------------------------------------
    call DiscretizedLapalacian(laplacian, N)
    call HarmonicPotential(harmonic_potential, N, dx, L, m, omega)
    H = -((hbar**2)/(2*m*dx**2))*laplacian + harmonic_potential
    ! --------------------------------------------------------------------------------------------

    ! ---- diagonalization -----------------------------------------------------------------------
    ! use the 'info' flag of the subroutine 'zheev' to check if the
    ! diagonalization has been successful
    info = 1
    do while (info .ne. 0)
        call ComputeEigenvalues(H, eig, info)
    end do 
    H = H / sqrt(dx) 
    ! --------------------------------------------------------------------------------------------

    ! ---- the starting state is the ground state ------------------------------------------------
    V_t(:, 1) = 0.d0
    psi_t(:, 1) = H(:, 1)
    prob_t(:, 1) = abs(psi_t(:, 1))**2 
    ! --------------------------------------------------------------------------------------------

    deallocate(eig, H, harmonic_potential, laplacian)

    ! ---- compute the time evolution of psi -----------------------------------------------------
    call psiTimeEvol(psi_t, prob_t, x_grid, p_grid, t_len, dt, omega, m, hbar, T, V_t)
    ! --------------------------------------------------------------------------------------------

    ! ---- compute the time evolution of the mean ------------------------------------------------
    do jj = 1, t_len 
        mean_t(jj) = sum(prob_t(:, 1) * prob_t(:, jj+1) * dx) 
    end do 
    print *, mean_t(9000)
    ! --------------------------------------------------------------------------------------------

    ! ---- save results --------------------------------------------------------------------------
    !call WriteRealMatrix(realpart(psi_t), x_grid, 'psi_real_time_evol_'//trim(str_r_d(T))//'.txt')
    !call WriteRealMatrix(imagpart(psi_t), x_grid, 'psi_imag_time_evol_'//trim(str_r_d(T))//'.txt')
    !call WriteRealMatrix(V_t, x_grid, 'V_time_evol_'//trim(str_r_d(T))//'.txt')
    call WriteRealMatrix(prob_t, x_grid, 'prob_time_evol_'//trim(str_r_d(T))//'.txt')
    call WriteRealVector(mean_t, 'pr_mean_time_evol_'//trim(str_r_d(T))//'.txt')
    ! --------------------------------------------------------------------------------------------
    
end program time_evolution