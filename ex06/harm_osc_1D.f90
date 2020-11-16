module HarmonicOscillator1D

    implicit none 
    
    ! variable to loop
    integer :: ii 

contains 

    subroutine DiscretizedLapalacian(lap, grid_size)

        ! this subroutine computes the discretized laplacian using
        ! the discretization psi'' = (psi_{n+1}+psi_{n-1}+2*psi_{n})/2
        
        ! laplacian matrix
        complex*16, dimension(:,:) :: lap 
        ! points in the grid 
        integer :: grid_size 

        ! the discretized laplacian has -2 in the diagonal
        ! and 1 in the sub-diagonals
        lap(1, 1) = -2
        do ii = 2, grid_size
            lap(ii, ii) = -2
            lap(ii-1, ii) = 1
            lap(ii, ii-1) = 1
        end do 

        return 

    end subroutine DiscretizedLapalacian


    subroutine HarmonicPotential(harmpot, N, dx, L, mass, omega)

        ! this subroutine computes the harmonic potential 

        ! matrix for the harmonic potential 
        complex*16, dimension(:,:) :: harmpot 
        ! omega, particle mass and discretization spacing 
        real*8 :: omega, mass, dx
        ! number of points 
        integer :: N, L 

        ! the harmonic potential is given by (m*omega^2*x^2)/2
        harmpot = 0.d0 
        do ii = 1, N 
            harmpot(ii, ii) = 0.5*mass*(omega*dx*(ii-L-1))**2
        end do 

        return

    end subroutine HarmonicPotential


    subroutine ComputeEigenvalues(matr, eig, info)

        ! this subroutine computes the eigenvalues of a complex hermitian matrix
        ! using the LAPACK subroutine 'zheev', which takes this arguments:
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


    subroutine ComputeProbabilityDens(H, prob) 

        ! this subroutine computes the probability densities 
        ! corresponding to the eigenvectors

        ! matrix which columns are the eigenvectors
        complex*16, dimension(:,:) :: H  
        ! matrix for the probability densities
        real*8, dimension(size(H, 1), size(H, 2)) :: prob

        ! pr = |psi|^2
        do ii = 1, size(H, 2)
            prob(:, ii) = abs(H(:, ii))**2
        end do 

        return 

    end subroutine ComputeProbabilityDens

end module HarmonicOscillator1D