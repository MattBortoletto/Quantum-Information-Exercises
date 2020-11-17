program time_evolution 

    use HarmonicOscillator1D
    use debugger 
    use Utilities 

    implicit none 

    complex*16, dimension(:), allocatable :: psi_0
    complex*16, dimension(:,:), allocatable :: psi_t 
    real*8, dimension(:,:), allocatable :: psi_RealPart, psi_ImPart
    integer :: L, t_len 
    real*8 :: dx, omega, m, hbar, x_min, x_max, t_min, t_max, T, dt 
    ! H: hamiltonian
    ! laplacian: discretized laplacian 
    ! harmonic_potential: harmonic potential
    complex*16, dimension(:,:), allocatable :: H, laplacian, harmonic_potential
    ! eig: eigenvalues array
    ! x_grid_points
    ! t_grid_points
    real*8, dimension(:), allocatable :: eig, x_grid_points, t_grid_points

    ! ---- initialize parameters ----
    L = 500
    t_len = 300  
    dx = 0.01
    t_min = 0.0
    t_max = 1.0 
    T = 3*t_max
    m = 1
    hbar = 1
    omega = 5 
    ! -------------------------------

    ! ---- compute the hamiltonian -------------------------------
    allocate(laplacian(N, N)) 
    allocate(harmonic_potential(N, N))
    call DiscretizedLapalacian(laplacian, N)
    call HarmonicPotential(harmonic_potential, N, dx, L, m, omega)

    H = -((hbar**2)/(2*m*dx**2))*laplacian + harmonic_potential
    ! ------------------------------------------------------------

    ! ---- diagonalization ----------------------------------------
    allocate(eig(N))
    ! use the 'info' flag of the subroutine 'zheev' to check if the
    ! diagonalization has been successful
    info = 1
    do while (info .ne. 0)
        call ComputeEigenvalues(H, eig, info)
    end do 
    ! -------------------------------------------------------------

    ! the starting state is the ground state 
    psi_0 = H(:, 1)
    
    deallocate(eig, H) ! harmonic_potential, laplacian 

    ! ---- compute the time evolution of psi ----
    allocate(psi_t(N, T_len))
    ! -------------------------------------------








end program time_evolution