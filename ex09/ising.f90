module Ising 

    implicit none 

    contains 

    subroutine IsingHamiltonian(H, N, lambda)  

        complex*16, dimension(:,:) :: H
        complex*16, dimension(:,:), allocatable :: offdiag  
        complex*16, dimension(:), allocatable :: diag
        integer :: N, ii, jj, kk, info 
        real*8 :: lambda 

        allocate(diag(2**N), offdiag(2**N, 2**N), stat=info)
        if (info .ne. 0) then 
            print *, "Allocation error!"
            stop 
        end if 

        do ii = 0, 2**N-1
            do jj = 0, N-1 
                diag(ii+1) = diag(ii+1) - 2*mod(ii/(2**jj), 2) 
            end do 
        end do
        
        do ii = 0, 2**N-1
            do jj = 0, 2**N-1
                offdiag(jj+1, ii+1) = 0
                do kk = 1, N-1 
                    if (ieor(ii, jj) == (2**(kk-1)+2**kk)) then 
                        offdiag(jj+1, ii+1) = offdiag(jj+1, ii+1) + 1
                    end if 
                end do 
            end do 
        end do 

        do ii = 1, 2**N 
            H(ii, ii) = lambda*diag(ii) + offdiag(ii, ii) 
        end do 

        deallocate(diag, offdiag)

        return 

    end subroutine IsingHamiltonian


    subroutine Diagonalize(matr, eig, info)

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

    end subroutine Diagonalize

end module Ising 