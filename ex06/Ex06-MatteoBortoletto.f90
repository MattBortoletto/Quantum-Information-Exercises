program harmonic_oscillator_1D

    use HarmonicOscillator1D 
    use Utilities
    use debugger 

    implicit none 

    ! N: number of points
    ! L: half of the interval L_s points
    ! info: 'zheev' subroutine flag
    ! jj: variable to loop
    integer :: N, L, info, jj
    ! dx: grid spacing
    ! omega: angular frequency of the oscillator
    ! m: mass
    ! hbar: Planck constant
    real*8 :: dx, omega, m, hbar
    ! H: hamiltonian
    ! laplacian: discretized laplacian 
    ! harmonic_potential: harmonic potential
    complex*16, dimension(:,:), allocatable :: H, laplacian, harmonic_potential
    ! eig: eigenvalues array
    ! grid_points: discretized grid points
    real*8, dimension(:), allocatable :: eig, grid_points
    ! probability densities array
    real*8, dimension(:,:), allocatable :: probabilities
    ! output file names
    character(:), allocatable :: energies_filename, states_filename, probabilities_filename
    ! flag to choose between default/custom parameters
    character(1) :: which_param
    ! flag to enable debugging 
    logical :: enable_debug 

    ! ---- enable/disable debugging ----
    enable_debug = .false.
    ! ----------------------------------
    
    ! ---- set the parameters ---------------------------------------------------------------------------------------
    print *, "Do you want to use use custom parameters of the default one (L=500, dx=0.01, omega=5, m=hbar=1)? [c/d]"
    read *, which_param
    if (which_param == "d") then 
        ! default values
        L = 500
        dx = 0.01
        omega = 5
        m = 1.0
        hbar = 1.0
    else if (which_param == "c") then 
        print *, "Please enter L, dx, omega, m and hbar:"
        read *, L, dx, omega, m, hbar 
    else 
        print *, "Invalid input."
        stop 
    end if 
    ! ----------------------------------------------------------------------------------------------------------------

    call checkpoint(debug=enable_debug, variable=L, message="L:", end_program=.false.)
    call checkpoint(debug=enable_debug, variable=dx, message="dx:", end_program=.false.)
    call checkpoint(debug=enable_debug, variable=omega, message="omega:", end_program=.false.)
    call checkpoint(debug=enable_debug, variable=m, message="m:", end_program=.false.)
    call checkpoint(debug=enable_debug, variable=hbar, message="hbar:", end_program=.false.)

    ! compute the number of points
    N = L*2 + 1

    call checkpoint(debug=enable_debug, variable=N, message="Number of points:", end_program=.false.)

    ! compute the values of the grid points 
    allocate(grid_points(N))
    do jj = 1, N 
        ! x_min = -L*dx
        grid_points(jj) = -L*dx + dx*(jj-1)
    end do 

    !call checkpoint(debug=(size(grid_points)==N), variable=N, message="grid_points dimension is correct!", end_program=.false.)

    ! ---- compute the hamiltonian -------------------------------
    allocate(laplacian(N, N)) 
    allocate(harmonic_potential(N, N))
    call DiscretizedLapalacian(laplacian, N)
    call HarmonicPotential(harmonic_potential, N, dx, L, m, omega)

    H = -((hbar**2)/(2*m*dx**2))*laplacian + harmonic_potential
    ! ------------------------------------------------------------

    allocate(eig(N))

    ! use the 'info' flag of the subroutine 'zheev' to check if the
    ! diagonalization has been successful
    info = 1
    do while (info .ne. 0)
        call ComputeEigenvalues(H, eig, info)
    end do 

    call checkpoint(debug=enable_debug, array_variable=eig, message="Eigenvalues:", end_program=.false.)

    ! normalize 
    H = H / sqrt(dx) 
    do ii = 1, N 
        H(:, ii) = H(:, ii) / norm2(real(H(:, ii)))
    end do

    call checkpoint(debug=enable_debug, matrix_variable=H, message="H:", end_program=.false.)

    ! compute the probability densities
    allocate(probabilities(N, N))
    call ComputeProbabilityDens(H=H, prob=probabilities)

    call checkpoint(debug=enable_debug, matrix_variable=probabilities, message="Probability densities:", end_program=.false.)

    ! save the results in text files with proper names
    energies_filename = "en_"//trim(str_i(L))//"_"//trim(str_r_e(dx))//"_"//trim(str_r_e(omega))//"_"&
                        &//trim(str_r_d(m))//"_"//trim(str_r_d(hbar))//".txt"

    states_filename = "ef_"//trim(str_i(L))//"_"//trim(str_r_e(dx))//"_"//trim(str_r_e(omega))//"_"&
                      &//trim(str_r_d(m))//"_"//trim(str_r_d(hbar))//".txt"

    probabilities_filename = "pr_"//trim(str_i(L))//"_"//trim(str_r_e(dx))//"_"//trim(str_r_e(omega))//"_"&
                             &//trim(str_r_d(m))//"_"//trim(str_r_d(hbar))//".txt"

    call WriteEigenvalues(eig, energies_filename)
    call WriteEigenvectors(real(H), grid_points, states_filename)
    call WriteEigenvectors(probabilities, grid_points, probabilities_filename)

    deallocate(H, grid_points, laplacian, harmonic_potential, eig, probabilities)

end program harmonic_oscillator_1D