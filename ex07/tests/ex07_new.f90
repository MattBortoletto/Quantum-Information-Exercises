subroutine Normalize(dx, N, psi, prob)

    integer :: N, ii 
    real*8 :: dx, normalization
    complex*16, dimension(N) :: psi 
    real*8, dimension(N) :: prob 

    normalization = sum(psi*conjg(psi))
    print*, 'NORMALIZATION: ',normalization
    do ii = 1, N 
        prob(ii) =  real(psi(ii)*conjg(psi(ii))/normalization)/dx 
    end do

end subroutine Normalize 

subroutine tev(dx, dp, dt, N, psi, x_max)

    use, intrinsic :: iso_c_binding

    implicit none 

    include 'fftw3.f03'

    real*8 :: dx, dp, dt, x_max, q0, T 
    integer :: N, ii 
    complex*16, dimension(N) :: psi, psi_fft
    integer*8 :: plan 

    T = 0.02 
    q0 = q0 - dt / T 

    print*, 'Q0 POSITION :', -q0

    do ii = 1, N
        psi(ii) = psi(ii)*exp(-0.5*dcmplx(0, (ii*dx-q0-x_max)**2)*dt) 
    end do

    ! Fourier transform
    call dfftw_plan_dft_1d(plan,N,psi,psi_fft,FFTW_FORWARD,FFTW_ESTIMATE)
    call dfftw_execute_dft(plan, psi, psi_fft)
    call dfftw_destroy_plan(plan)    

    ! step : Kinetic
    do ii = 1, N/2
        psi_fft(ii) = psi_fft(ii)*exp(-0.5*dcmplx(0, ((ii-1)*dp)**2)*dt) 
    end do

    do ii = N/2+1, N
        psi_fft(ii) = psi_fft(ii)*exp(-0.5*dcmplx(0, ((ii-1-N)*dp)**2)*dt) 
    end do


    ! Fourier antitransform
    call dfftw_plan_dft_1d(plan,N,psi_fft,psi,FFTW_BACKWARD,FFTW_ESTIMATE)
    call dfftw_execute_dft(plan, psi_fft, psi)
    call dfftw_destroy_plan(plan)  

    ! xmax step : Potential
    do ii = 1, N
        psi(ii) = psi(ii)*exp(-0.5*dcmplx(0, (ii*dx-q0-x_max)**2)*dt)/N 
    end do

end subroutine tev 



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

    ! ---- normalization + probability -------------
    call Normalize(dx, N, psi_t(:, 1), prob_t(:, 1))
    ! ----------------------------------------------

    ! ---- compute the time evolution of psi ------------------------------------
    tmp_ps = psi_t(:, 1)
    tmp_pr = prob_t(:, 1)
    do jj = 1, NT 
        call tev(dx, dp, dt, N, tmp_ps, x_max)
        call Normalize(dx, N, tmp_ps, tmp_pr)
        psi_t(:, jj+1) = tmp_ps 
        prob_t(:, jj+1) = tmp_pr 
    end do 
    ! ---------------------------------------------------------------------------

    ! ---- grid -------------------------
    do jj = 1, N 
        x_grid(jj) = - x_max + (jj-1)*dx 
    end do 
    ! ----------------------------------

    ! ---- save results ---------------------------------------------------
    ! call WriteRealMatrix(realpart(psi_t), x_grid, 'psi_real_time_evol.txt')
    ! call WriteRealMatrix(imagpart(psi_t), x_grid, 'psi_imag_time_evol.txt')
    call WriteRealMatrix(prob_t, x_grid, 'prob_time_evol.txt')
    ! ---------------------------------------------------------------------

    print *, "END!"

    deallocate(eig, H, harmonic_potential, laplacian)
    
end program time_evolution