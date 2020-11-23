program time_evolution 

    use HarmonicOscillator1D
    use debugger 
    use Utilities 
    use TimeEvolution 

    implicit none 

    ! complex matrix for the time evolution of the psi
    ! the first column contains the x grid points and the
    ! other columns the psi at the different time steps
    complex*16, dimension(:,:), allocatable :: psi_t
    ! L: semi-length of the space interval (number of points)
    ! t_len: number of points for t
    ! N: number of points for x 
    ! info: stats flag
    ! jj: variable to loop
    integer :: L, t_len, N, info, jj 
    ! dx: discretization for x
    ! omega: angular frequency
    ! m: mass 
    ! hbar: Planck constant
    ! T: hamiltonian constant
    ! dt: discretization for t
    ! pi: greek pi
    real*8 :: dx, omega, m, hbar, T, dt, pi 
    ! H: hamiltonian
    ! laplacian: discretized laplacian 
    ! harmonic_potential: harmonic potential
    complex*16, dimension(:,:), allocatable :: H, laplacian, harmonic_potential
    ! eig: eigenvalues array
    ! x_grid: x points
    ! p_grid: p points
    ! mean_t: mean position of the pdf at time t
    real*8, dimension(:), allocatable :: eig, x_grid, p_grid, mean_t
    ! prob_t: probability density at time t. The first column contains 
    !         the x grid and the others the time evolution of the probability 
    !         density
    ! V_t: time evolution of the potential. The first column contains the x
    !      grid and the others the time evolution of the potential
    real*8, dimension(:,:), allocatable :: prob_t, V_t 

    ! ---- initialize parameters -----------------------------------------------------------------
    L = 500
    dx = 0.01
    N = 2*L + 1

    print *, "Please enter the value of T: "
    read *, T

    dt = 0.01
    t_len = int( T / dt )
    print *, t_len 
    m = 1
    hbar = 1
    omega = 1
    pi = 4.d0 * datan(1.d0)
    ! --------------------------------------------------------------------------------------------

    ! ---- allocate memory -----------------------------------------------------------------------
    allocate(x_grid(N), p_grid(N), laplacian(N, N), harmonic_potential(N, N), eig(N),& 
             &psi_t(N, t_len+1), prob_t(N, t_len+1), V_t(N, t_len+1), mean_t(t_len), stat=info)
    call Checkpoint(debug=info .ne. 0, message="Memory allocation fail!", end_program=.true.)
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
    !call Checkpoint(debug=.false., array_variable=x_grid, message="x_grid:", end_program=.false.)
    !call Checkpoint(debug=.false., array_variable=p_grid, message="p_grid:", end_program=.false.)
    ! --------------------------------------------------------------------------------------------

    ! ---- compute the hamiltonian ---------------------------------------------------------------
    call DiscretizedLapalacian(laplacian, N)
    call HarmonicPotential(harmonic_potential, N, dx, L, m, omega)
    H = -((hbar**2)/(2*m*dx**2))*laplacian + harmonic_potential
    ! --------------------------------------------------------------------------------------------

    ! ---- diagonalization -----------------------------------------------------------------------
    ! use the 'info' flag of the subroutine 'zheev' to check if the
    ! diagonalization has been successful
    call ComputeEigenvalues(H, eig, info)
    call Checkpoint(debug=info .ne. 0, message="Diagonalization failed!", end_program=.true.)
    H = H / sqrt(dx) 
    ! --------------------------------------------------------------------------------------------

    ! ---- the starting state is the ground state ------------------------------------------------
    V_t(:, 1) = 0.d0
    psi_t(:, 1) = H(:, 1)
    prob_t(:, 1) = abs(psi_t(:, 1))**2 
    ! --------------------------------------------------------------------------------------------

    deallocate(eig, H, harmonic_potential, laplacian)

    ! ---- compute the time evolution of psi -----------------------------------------------------
    call psiTimeEvol(psi_t, prob_t, x_grid, p_grid, t_len, dt, omega, m, hbar, T, V_t, dx)
    ! --------------------------------------------------------------------------------------------

    ! ---- compute the time evolution of the mean ------------------------------------------------
    mean_t = 0.d0 
    do jj = 1, t_len 
        mean_t(jj) = sum(prob_t(:, jj) * x_grid * dx)
    end do 
    ! --------------------------------------------------------------------------------------------

    ! ---- save results --------------------------------------------------------------------------
    call WriteRealMatrix(realpart(psi_t), x_grid, 'psi_real_time_evol_'//trim(str_r_d(T))//'.txt')
    call WriteRealMatrix(imagpart(psi_t), x_grid, 'psi_imag_time_evol_'//trim(str_r_d(T))//'.txt')
    call WriteRealMatrix(V_t, x_grid, 'V_time_evol_'//trim(str_r_d(T))//'.txt')
    call WriteRealMatrix(prob_t, x_grid, 'prob_time_evol_'//trim(str_r_d(T))//'.txt')
    call WriteRealVector(mean_t, 'pr_mean_time_evol_'//trim(str_r_d(T))//'.txt')
    ! --------------------------------------------------------------------------------------------

    deallocate(psi_t, prob_t, x_grid, p_grid, mean_t, V_t)
    
end program time_evolution