subroutine shift_pot(Tmax, t, x_min, dt, dx, hbar, omega, m)

    real*8 :: hbar, omega, m, dx, dt, x_min, t, t_max
    
    fact = complex(0, -dt*0.25*m*omega**2 / hbar)
    x_sh = x_min - t/Tmax
    do ii = 1, N 
        psi(ii) = psi(ii) * exp(-fact*x_sh**2)
        x_sh = x_sh + dx 
    end do 

end subroutine shift_pot

subroutine tep()

    temp = psi0
    t = 0.0 
    do idx = 2, NT 
        
end subroutine tep 


program time_evolution 

    use HarmonicOscillator1D
    use debugger 
    use Utilities 
    use TimeEvolution 

    implicit none 

    complex*16, dimension(:,:), allocatable :: psi_t
    complex*16, dimension(:), allocatable :: tmp_ps
    integer :: NT, N, info, jj, L 
    real*8 :: dx, omega, m, hbar, dt, pi, x_max, dp
    complex*16, dimension(:,:), allocatable :: H, laplacian, harmonic_potential
    real*8, dimension(:), allocatable :: eig, tmp_pr, x_grid 
    real*8, dimension(:,:), allocatable :: prob_t

    print *, "NEW PROGRAM"

    ! ---- initialize parameters ----
    L = 500
    N = 2*L + 1 
    NT = 500
    dx = 0.01
    dt = 0.002
    x_max = dx*N/2
    dp = pi/x_max 

    m = 1
    hbar = 1
    omega = 1
    pi = 4.d0 * datan(1.d0)
    ! -------------------------------

    ! ---- allocate memory -----------
    allocate(laplacian(N, N)) 
    allocate(harmonic_potential(N, N))
    allocate(eig(N))
    allocate(psi_t(N, NT+1))
    allocate(prob_t(N, NT+1))
    allocate(tmp_ps(N))
    allocate(tmp_pr(N))
    allocate(x_grid(N))
    ! --------------------------------

    ! ---- compute the hamiltonian -------------------------------
    call DiscretizedLapalacian(laplacian, N)
    call HarmonicPotential(harmonic_potential, N, dx, L, m, omega)
    H = -((hbar**2)/(2*m*dx**2))*laplacian + harmonic_potential
    ! ------------------------------------------------------------

    ! ---- diagonalization ----------------------------------------
    info = 1
    do while (info .ne. 0)
        call ComputeEigenvalues(H, eig, info)
    end do 
    ! -------------------------------------------------------------

    ! ---- starting state 
    psi_t(:, 1) = H(:, 1)
    ! -------------------

    ! ---- compute the time evolution of psi ------------------------------------

    ! ---------------------------------------------------------------------------

 

    ! ---- save results ---------------------------------------------------
    ! call WriteRealMatrix(realpart(psi_t), x_grid, 'psi_real_time_evol.txt')
    ! call WriteRealMatrix(imagpart(psi_t), x_grid, 'psi_imag_time_evol.txt')
    call WriteRealMatrix(prob_t, x_grid, 'prob_time_evol.txt')
    ! ---------------------------------------------------------------------

    print *, "END!"

    deallocate(eig, H, harmonic_potential, laplacian)
    
end program time_evolution