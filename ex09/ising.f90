module Ising 

    implicit none 

    contains 

    subroutine IsingHamiltonian(H, N, lambda)  

        ! computes the Hamiltonian for the one-dimensional quantum 
        ! Ising model 

        ! Hamiltonian
        complex*16, dimension(:,:) :: H
        ! diagonal terms
        complex*16, dimension(:), allocatable :: diag
        ! N: number of subsystems
        ! ii, jj: variables to loop 
        ! info: stat flag
        integer :: N, ii, jj, info  
        ! lambda
        real*8 :: lambda 

        allocate(diag(2**N), stat=info)
        if (info .ne. 0) then 
            print *, "Allocation error!"
            stop 
        end if 

        ! first fill the diagonal terms
        diag = N 
        do ii = 0, 2**N-1
            do jj = 0, N-1 
                diag(ii+1) = diag(ii+1) - 2*mod(ii/(2**jj), 2)
            end do 
            diag(ii+1) = diag(ii+1) * lambda 
        end do
        
        ! then compute the interaction terms 
        do ii = 0, 2**N-1
            do jj = 0, 2**N-1
                H(jj+1, ii+1) = ComputeInter(jj, ii, N) 
            end do 
        end do 

        ! add the two terms to find the complete Hamiltonian
        do ii = 1, 2**N 
            H(ii, ii) = diag(ii) + H(ii, ii) 
        end do 

        deallocate(diag)

        return 

    end subroutine IsingHamiltonian


    function ComputeInter(q, p, N) result(res) 

        ! computes the interaction term, which 
        ! acts like a XOR, at position (q,p)

        integer :: q, p, N, ii, res 

        res = 0
        do ii = 1, N-1
            if ((2**(ii-1)+2**ii) == xor(q,p)) then 
                res = res + 1
            end if 
        end do 

        return 

    end function ComputeInter 


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