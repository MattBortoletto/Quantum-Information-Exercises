program time_evolution 

    use HarmonicOscillator1D
    use debugger 
    use Utilities 
    use TimeEvolution 

    implicit none 

    complex*16, dimension(:,:), allocatable :: psi_t
    integer :: L, t_len, N, info, jj 
    real*8 :: dx, omega, m, hbar, T, dt, pi 
    ! H: hamiltonian
    ! laplacian: discretized laplacian 
    ! harmonic_potential: harmonic potential
    complex*16, dimension(:,:), allocatable :: H, laplacian, harmonic_potential
    ! eig: eigenvalues array
    ! x_grid
    ! p_grid
    real*8, dimension(:), allocatable :: eig, x_grid, p_grid 


    ! ---- initialize parameters ----
    L = 500
    t_len = 300  
    dx = 0.01
    N = 2*L + 1
    T = 3.0
    dt = T / t_len ! t is in [0:T]
    m = 1
    hbar = 1
    omega = 5
    pi = 4.d0 * datan(1.d0)
    ! -------------------------------

    allocate(x_grid(N))
    allocate(p_grid(N))
    allocate(laplacian(N, N)) 
    allocate(harmonic_potential(N, N))
    allocate(eig(N))
    allocate(psi_t(N, t_len+1))

    ! ---- compute the grid points ----
    ! x grid  
    do jj = 1, N 
        ! x_min = -L*dx
        x_grid(jj) = -L*dx + dx*(jj-1)
    end do
    ! p grid 
    do jj = 1, N 
        ! lower = 0 
        ! upper = (2*pi*N)/(dx*(N-1))
        ! dp = (upper - lower) / (N-1) = upper / (N-1)
        p_grid(jj) = (jj-1)*(2*pi) / (dx*N)  
        p_grid(N/2:) = p_grid(N/2:) - p_grid(N) - p_grid(2)
    end do 
    
    ! print *, p_grid 

    ! ---- compute the hamiltonian -------------------------------
    call DiscretizedLapalacian(laplacian, N)
    call HarmonicPotential(harmonic_potential, N, dx, L, m, omega)

    H = -((hbar**2)/(2*m*dx**2))*laplacian + harmonic_potential
    ! ------------------------------------------------------------

    ! ---- diagonalization ----------------------------------------
    ! use the 'info' flag of the subroutine 'zheev' to check if the
    ! diagonalization has been successful
    info = 1
    do while (info .ne. 0)
        call ComputeEigenvalues(H, eig, info)
    end do 
    H = H / sqrt(dx) 
    ! -------------------------------------------------------------

    ! the starting state is the ground state 
    psi_t(:, 1) = H(:, 1)
    
    deallocate(eig, H) ! harmonic_potential, laplacian 

    ! ---- compute the time evolution of psi -------------------
    call psiTimeEvol(psi_t, x_grid, p_grid, t_len, dt, omega, m, hbar)
    ! ----------------------------------------------------------

    call WriteRealMatrix(realpart(psi_t), x_grid, 'psi_real_time_evol.txt')
    call WriteRealMatrix(imagpart(psi_t), x_grid, 'psi_imag_time_evol.txt')

    print *, "END!"
    
end program time_evolution