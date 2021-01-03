module Ising 

    use debugger

    implicit none 

    contains 

    subroutine IsingHamiltonian(H, N, lambda)  

        ! computes the Hamiltonian for the one-dimensional quantum 
        ! Ising model 

        ! Hamiltonian
        real*8, dimension(:,:) :: H
        ! diagonal terms
        real*8, dimension(:), allocatable :: diag
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


    subroutine DiagonalizeR(matr, eig, info)

        ! this subroutine computes the eigenvalues of a real symmetric matrix
        ! using the LAPACK subroutine 'dsyev', which takes this arguments:
        ! - jobz = "V": eigenvalues and eigenvectors are computed
        ! - uplo = "U": upper triangle of A is stored
        ! - N = size(matr, 1): order of the matrix
        ! - a = matr: input matrix
        ! - lda = size(matr, 1): leading dimension of the matrix
        ! - w: if info = 0, the eigenvalues in ascending order
        ! - work: complex array of dimension (max(1,lwork))
        !         on exit, if info = 0, work(1) returns the optimal lwork
        ! - lwork: length of the array work
        ! - info: output, if info = 0 the exit is successful

        real*8, dimension(:, :), intent(in) :: matr
        real*8, dimension(size(matr, 1)) :: eig 
        complex*16, dimension(:), allocatable :: work
        character(1) :: jobz, uplo
        integer :: n, lda, lwork, info

        n = size(matr, 1)
        lda = size(matr, 1)
        jobz = "V"
        uplo = "U"

        ! ------------ find the optimal lwork ------------------------------
        lwork = -1
        allocate(work(1))
        call dsyev(jobz, uplo, n, matr, lda, eig, work, lwork, info)
        ! on exit, if info = 0, work(1) returns the optimal lwork
        lwork = int(work(1))
        deallocate(work) 
        ! ------------------------------------------------------------------

        ! allocate work using the optimal lwork 
        allocate(work(lwork))
        
        ! perform the diagonalization
        call dsyev(jobz, uplo, n, matr, lda, eig, work, lwork, info)

        deallocate(work)

        return 

    end subroutine DiagonalizeR


    function KroneckerProd(m1, m2) result(kprod)

        ! computes the tensor product of two matrices (kprod = m1 x m2)

        real*8, dimension(:,:) :: m1, m2
        real*8, dimension(size(m1,1)*size(m2,1), size(m1,2)*size(m2,2)) :: kprod
        integer :: ii, jj 

        kprod = 0.d0 
        forall (ii = 1:size(m1,1), jj = 1:size(m1,2)) 
            kprod(size(m2,1)*(ii-1)+1 : size(m2,1)*ii , size(m2,2)*(jj-1)+1 : size(m2,2)*jj) = m1(ii,jj)*m2
        end forall 

        return 

    end function KroneckerProd 


    function idmatr(dim) result(id) 

        ! builds a (dim x dim) identity matrix

        integer :: ii, dim 
        real*8, dimension(dim, dim) :: id 

        id = 0.d0
        do ii = 1, dim
            id(ii, ii) = 1.d0
        end do 

        return 

    end function idmatr


    function RSRG(H_N, N, niter) result(gs)

        ! Real space Renormalization Group 

        real*8, dimension(:,:) :: H_N
        real*8, dimension(:,:), allocatable :: H_2N, H_2N_L, H_2N_R, Htmp
        integer :: N, niter, ii, info 
        real*8, dimension(:,:), allocatable :: P
        real*8, dimension(2,2) :: sigma_x
        real*8, dimension(:), allocatable :: eig
        real*8 :: gs 

        sigma_x(1, 1) = 0
        sigma_x(1, 2) = 1
        sigma_x(2, 1) = 1
        sigma_x(2, 2) = 0

        allocate(H_2N(2**(2*N), 2**(2*N)), &
                 H_2N_L(2**N, 2**N), &
                 H_2N_R(2**N, 2**N), &
                 Htmp(2**(2*N), 2**(2*N)), &
                 eig(2**(2*N)), stat=info)
        call checkpoint(debug=(info.ne.0), message="Allocation failed!", &
                        end_program=.true.)

        ! build H_2N
        H_2N_L = KroneckerProd(idmatr(2**(N-1)), sigma_x)
        H_2N_R = KroneckerProd(sigma_x, idmatr(2**(N-1)))
        H_2N = KroneckerProd(H_N,  idmatr(2**N)) + KroneckerProd(idmatr(2**N), H_N) + KroneckerProd(H_2N_L, H_2N_R)

        do ii = 1, niter 

            Htmp = H_2N 

            ! diagonalize 
            call DiagonalizeR(H_2N, eig, info)
            call checkpoint(debug=(info.ne.0), message="Diagonalization failed!", &
                            end_program=.true.)
                            
            ! build the projector
            P = H_2N(:, 1:2**N) 

            H_2N = Htmp 

            ! project: Htilde_N = P^+ H_2N P
            H_N = matmul(transpose(P), matmul(H_2N, P))
            
            ! build Htilde_2N = HtildeN_L + HtildeN_R + HtildeN_int
            !                   HtildeN*1 + 1*HtildeN + Atilde_N*Btilde_N
            ! where Atilde_N = P^+ A P, Btilde_N = P^+ B P 
            H_2N_L = matmul(transpose(P), matmul(KroneckerProd(idmatr(2**N), H_2N_L), P))
            H_2N_R = matmul(transpose(P), matmul(KroneckerProd(H_2N_R, idmatr(2**N)), P))

            H_2N = KroneckerProd(H_N, idmatr(2**N)) + KroneckerProd(idmatr(2**N), H_N) + KroneckerProd(H_2N_L, H_2N_R)
            
        end do 

        ! diagonalize 
        call DiagonalizeR(H_2N, eig, info)
        call checkpoint(debug=(info.ne.0), message="Diagonalization failed!", &
                        end_program=.true.)
        
        gs = eig(1)

        deallocate(H_2N, H_2N_L, H_2N_R, Htmp, P, eig)

        return 

    end function RSRG 

end module Ising 