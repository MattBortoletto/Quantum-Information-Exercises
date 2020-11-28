module DensityMatrices

    implicit none 

contains 

    !subroutine InitPureSepState(D, N, state)
    
    function ComplexSquareMatrixTrace(matrix) result(trace) 

        ! computes the trace of a square complex matrix 

        complex*16, dimension(:,:) :: matrix 
        complex*16 :: trace 
        integer :: ii 

        if (size(matrix, 1) .ne. size(matrix, 2)) then 
            print *, "Error! Non square matrix."
        else 
            trace = 0.d0
            do ii = 1, size(matrix, 1)
                trace = trace + matrix(ii, ii)
            end do 
        end if 

    end function ComplexSquareMatrixTrace


    subroutine BuildSeparableStateFromMatr(coeff_matrix, psi) 

        ! constructs a separable state starting from the coefficients matrix

        integer :: ii, jj, D, N 
        complex*16, dimension(:,:) :: coeff_matrix 
        complex*16, dimension(size(coeff_matrix, 1)**size(coeff_matrix, 2)) :: psi 

        D = size(coeff_matrix, 1)
        N = size(coeff_matrix, 2)

        psi = dcmplx(1d0) 
        do ii = 1, D**N 
            do jj = 1, N 
                print *, mod( (ii-1) / (D**(N-jj)), D)
                psi(ii) = psi(ii) * coeff_matrix(mod( (ii-1) / (D**(N-jj)) , D) + 1, jj)
            end do
        end do 

    end subroutine BuildSeparableStateFromMatr


    subroutine BuildSeparableStateFromVect(coeff_vect, psi, N, D)  

        ! constructs a separable state starting from the coefficients vector

        integer :: ii, jj, N, D 
        complex*16, dimension(:) :: coeff_vect 
        complex*16, dimension(:) :: psi 

        do ii = 1, D**N 
            psi(ii) = dcmplx(1d0) 
            do jj = 1, N 
                print *, mod( (ii-1) / (D**(N-jj)), D)
                psi(ii) = psi(ii) * coeff_vect(mod( (ii-1) / (D**(N-jj)) , D) + 1)
            end do
        end do 

    end subroutine BuildSeparableStateFromVect
            
    
    subroutine ReducedDensityMatrixB(rho, D, rho_B) 

        ! computes the reduced density matrix for N=2

        integer :: D, jj, kk, ll
        complex*16, dimension(:,:) :: rho 
        ! the reduced density matrix has dimension D x D
        complex*16, dimension(D, D) :: rho_B

        rho_B = 0.d0
        do jj = 1, D  
            do kk = jj, D 
                do ll = 1, D 
                    rho_B(jj, kk) = rho_B(jj, kk) + rho((ll-1)*D+jj, (ll-1)*D+kk)
                end do 
                if (jj /= kk) rho_B(kk, jj) = conjg(rho_B(jj, kk))
            end do 
        end do 

    end subroutine ReducedDensityMatrixB


    subroutine ReducedDensityMatrixA(rho, D, rho_A) 

        ! computes the reduced density matrix for N=2

        integer :: D, jj, kk, ll
        complex*16, dimension(:,:) :: rho 
        ! the reduced density matrix has dimension D x D
        complex*16, dimension(D, D) :: rho_A

        rho_A = 0.d0
        do jj = 1, D  
            do kk = jj, D 
                do ll = 1, D
                    rho_A(jj, kk) = rho_A(jj, kk) + rho((jj-1)*D+ll, (kk-1)*D+ll)
                end do 
                if (jj /= kk) rho_A(kk, jj) = conjg(rho_A(jj, kk))
            end do 
        end do 

    end subroutine ReducedDensityMatrixA


    subroutine WriteComplex16Matr(matr, trace, filename) 

        complex*16, dimension(:,:) :: matr
        complex*16 :: trace 
        character(*) :: filename 
        integer :: ll

        open(unit = 73, file = filename, action = "write", status = "replace") 
        
        write(73, *) "MATRIX:"
        do ll = 1, size(matr, 1)
            write(73, *) matr(ll, :) 
        end do
        write(73, *) " "
        ! write(73, *) "DIMENSION: rows =",  M%dim(1), ",  columns = ", M%dim(2)
        ! write(73, *) " "
        write(73, *) "TRACE:", trace 
        write(73, *) " "
        close(73)

        return 

    end subroutine WriteComplex16Matr

end module DensityMatrices