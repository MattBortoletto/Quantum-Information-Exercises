module HarmonicOscillator1D

    implicit none 
    
    integer :: ii 

contains 

    subroutine DiscretizedLapalacian(lap, grid_size, periodic_bc)
        
        complex*16, dimension(:,:) :: lap 
        integer :: grid_size 
        logical :: periodic_bc 

        lap(1, 1) = -2

        do ii = 2, grid_size
            lap(ii, ii) = -2
            lap(ii-1, ii) = 1
            lap(ii, ii-1) = 1
        end do 

        if (periodic_bc) then 
            lap(1, grid_size) = 1
            lap(grid_size, 1) = 1
        end if 

        return 

    end subroutine DiscretizedLapalacian


    subroutine HarmonicPotential(harmpot, N, dx, L, mass, omega)

        complex*16, dimension(:,:) :: harmpot 
        real*8 :: omega, mass, dx
        integer :: N, L 

        harmpot = 0.d0 
        do ii = 1, N 
            harmpot(ii, ii) = 0.5*mass*(omega*dx*(ii-L-1)) 
        end do 

        return

    end subroutine HarmonicPotential


    subroutine ComputeEigenvalues(matr, eig, info)

        ! this subroutine computes the eigenvalues of a complex hermitian matrix
        ! using the LAPACK subroutine 'cheev', which takes this arguments:
        ! - jobz = "V": eigenvalues and eigenvectors are computed
        ! - uplo = "U": upper triangle of A is stored
        ! - N = size(matr, 1): order of the matrix
        ! - a = matr: input matrix
        ! - lda = size(matr, 1): leading dimension of the matrix
        ! - w: if info = 0, the eigenvalues in ascending order
        ! - work: complex array of dimension (max(1,lwork))
        !         on exit, if info = 0, work(1) returns the optimal lwork
        ! - lwork: length of the array work
        ! - rwork: workspace
        !          real array, dimension (max(1, 3*n-2))
        ! - info: output, if info = 0 the exit is successful

        complex*16, dimension(:, :), intent(in) :: matr
        real*8, dimension(size(matr, 1)) :: eig 
        real*8, dimension(:), allocatable :: rwork
        complex*16, dimension(:), allocatable :: work
        character(1) :: jobz, uplo
        integer :: n, lda, lwork, info

        n = size(matr, 1)
        lda = size(matr, 1)
        jobz = "V"
        uplo = "U"
        allocate(rwork(max(1, 3*size(matr, 1)-2)))   

        ! ------------ find the optimal lwork ------------------------------
        lwork = -1
        allocate(work(1))
        call zheev(jobz, uplo, n, matr, lda, eig, work, lwork, rwork, info)
        ! on exit, if info = 0, work(1) returns the optimal lwork
        lwork = int(work(1))
        deallocate(work) 
        ! ------------------------------------------------------------------

        ! allocate work using the optimal lwork 
        allocate(work(lwork))
        
        ! perform the diagonalization
        call zheev(jobz, uplo, n, matr, lda, eig, work, lwork, rwork, info)

        deallocate(rwork, work)

        return 

    end subroutine ComputeEigenvalues


    subroutine ComputeProb(H, prob) 

        complex*16, dimension(:,:) :: H  
        real*8, dimension(size(H, 1), size(H, 2)) :: prob

        do ii = 1, size(H, 2)
            prob(:, ii) = abs(H(:, ii))**2
        end do 

        return 

    end subroutine ComputeProb

end module HarmonicOscillator1D


module Utilities

    implicit none 

contains 

    function str(k) result(k_str)

        ! this function is used to convert an integer to string

        class(*), intent(in) :: k
        character(5) :: k_str

        select type(k) 
            type is (integer(2))
                write (k_str, *) k
                k_str = adjustl(k_str)
            type is (integer(4))
                write (k_str, *) k
                k_str = adjustl(k_str)
            type is (real(4))
                write (k_str, *) k
                k_str = adjustl(k_str)
            type is (real(8))
                write (k_str, *) k
                k_str = adjustl(k_str)
        end select

        return 

    end function str


    subroutine WriteEigenvalues(eig, filename) 

        character(*) :: filename
        real*8, dimension(:) :: eig 
        integer :: ii 

        open(unit=73, file=filename, action="write", status="replace")
        do ii = 1, size(eig) 
            write(73, *) eig(ii)
        end do 
        close(73)

        return 

    end subroutine WriteEigenvalues


    subroutine WriteEigenvectors(H, prob, grid_points, filename) 

        character(*) :: filename
        real*8, dimension(:,:) :: H, prob  
        real*8, dimension(:) :: grid_points  
        integer :: ii, jj 

        open(unit=73, file=filename, action="write", status="replace")
        do ii = 1, size(prob, 1) 
            do jj = 1, size(prob, 2)
                write(73, *) grid_points(ii), H(ii, jj), prob(ii, jj)
            end do 
        end do 
        close(73)

        return 
        
    end subroutine WriteEigenvectors

end module Utilities



program harmonic_oscillator_1D

    use HarmonicOscillator1D 
    use Utilities

    implicit none 

    integer :: N, L, info, jj
    real*8 :: dx, omega, m, hbar
    complex*16, dimension(:,:), allocatable :: H, laplacian, harmonic_potential
    real*8, dimension(:), allocatable :: eig, grid_points
    real*8, dimension(:,:), allocatable :: probabilities
    character(:), allocatable :: energies_filename, states_filename
    
    ! ---- default values ----
    L = 500
    dx = 1.0e-04
    omega = 1.0e04
    m = 1.0
    hbar = 1.0
    ! ------------------------

    N = L*2 + 1

    allocate(grid_points(N))
    do jj = 1, N 
        ! x_min = -L*dx
        grid_points(jj) = -L*dx + dx*(jj-1)
    end do 

    allocate(laplacian(N, N)) 
    allocate(harmonic_potential(N, N))
    call DiscretizedLapalacian(laplacian, N, .false.)
    call HarmonicPotential(harmonic_potential, N, dx, L, m, omega)

    H = -((hbar**2)/(2*m*dx**2))*laplacian + harmonic_potential

    allocate(eig(N))
    ! use the info flag of the subroutine 'cheev' to check if the
    ! diagonalization has been successful
    info = 1
    do while (info .ne. 0)
        call ComputeEigenvalues(H, eig, info)
    end do 

    ! normalize 
    H = H / sqrt(dx)

    allocate(probabilities(N, N))
    call ComputeProb(H, probabilities)

    energies_filename = "e_"//trim(str(L))//"_"//trim(str(dx))//"_"//trim(str(omega))//"_"&
                        &//trim(str(m))//"_"//trim(str(hbar))//".txt"

    states_filename = "s_"//trim(str(L))//"_"//trim(str(dx))//"_"//trim(str(omega))//"_"&
                      &//trim(str(m))//"_"//trim(str(hbar))//".txt"

    call WriteEigenvalues(eig, energies_filename)
    call WriteEigenvectors(real(H), probabilities, grid_points, states_filename)

    deallocate(grid_points, laplacian, harmonic_potential, eig, probabilities)

end program harmonic_oscillator_1D