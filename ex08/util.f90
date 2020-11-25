module Utilities

    use debugger

    implicit none 

contains 

    function str_i(k) result(str)

        ! this function converts an integer into string

        character(len=20) :: str
        integer, intent(in) :: k

        write (str, *) k
        str = adjustl(str)

        return 

    end function str_i


    function str_r_d(k) result(str)

        ! this function converts a real into string 
        ! (decimal notation)

        character(len=20) :: str
        real*8, intent(in) :: k

        write (str, "(F10.5)") k
        str = adjustl(str)

        return 

    end function str_r_d


    function str_r_e(k) result(str)
        
        ! this function converts a real into string 
        ! (exponential notation)

        character(len=20) :: str 
        real*8, intent(in) :: k

        write (str, "(E10.1)") k
        str = adjustl(str)

        return 

    end function str_r_e


    subroutine WriteRealVector(vect, filename) 

        ! this subroutine writes a vector into an output file

        character(*) :: filename
        real*8, dimension(:) :: vect 
        integer :: ii 

        open(unit=73, file=filename, action="write", status="replace")
        do ii = 1, size(vect) 
            write(73, *) ii-1, vect(ii)
        end do 
        close(73)

        return 

    end subroutine WriteRealVector


    subroutine WriteRealMatrix(matr, grid_points, filename) 

        ! this subroutine writes a matrix into an output file

        character(*) :: filename
        real*8, dimension(:,:) :: matr 
        real*8, dimension(:) :: grid_points  
        integer :: ii

        open(unit=73, file=filename, action="write", status="replace")
        do ii = 1, size(matr, 1) 
                write(73, *) grid_points(ii), matr(ii, :) 
        end do 
        close(73)

        return 
        
    end subroutine WriteRealMatrix


    subroutine NormalizePsi(psi, dx)

        ! this subroutine normalizes the psi and then checks that
        ! the normalization has been done properly

        ! the Checkpoint subroutine is contained in the 'debugger' module
    
        complex*16, dimension(:), allocatable :: psi
        real*8 :: dx, norm
        integer :: ii
        
        ! normalize 
        norm = 0.0
        do ii = 1, size(psi)
           norm = norm + (real(psi(ii))**2 + aimag(psi(ii))**2) * dx
        end do
    
        psi = psi / sqrt(norm)
    
        ! check 
        norm = 0.0
        do ii = 1, size(psi)
            norm = norm + (real(psi(ii))**2 + aimag(psi(ii))**2) * dx
        end do
    
        if (abs(norm-1.0) .gt. 0.0001) then
            print *, "Normalization error!"
            stop
        end if  
    
    end subroutine NormalizePsi


    subroutine OuterProd(a, b, aob) 
        
        ! this subroutine computes the outer product of two vectors

        complex*16, dimension(:) :: a, b 
        complex*16, dimension(:,:) :: aob 

        aob = spread(a, dim=2, ncopies=size(b)) * spread(b, dim=1, ncopies=size(a)) 

    end subroutine OuterProd

    
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


    subroutine BuildSeparableState(coeff_matrix, psi) 

        ! constructs a separable state starting from the coefficients matrix

        integer :: ii, jj 
        complex*16, dimension(:,:) :: coeff_matrix 
        complex*16, dimension(size(coeff_matrix, 1)**size(coeff_matrix, 2)) :: psi 
        
        integer :: D, N 
        integer, dimension(size(coeff_matrix, 2)) :: den 

        D = size(coeff_matrix, 1)
        N = size(coeff_matrix, 2)
        
        psi = dcmplx(1d0)
        ! do ii = 1, size(coeff_matrix, 1)**size(coeff_matrix, 2) ! 1, D**N 
        !     do jj = 1, size(coeff_matrix, 2) ! 1, N 
        !         psi(ii) = psi(ii) * &
        !                   coeff_matrix(mod((ii-1)/(int(size(coeff_matrix, 1)**(size(coeff_matrix, 2)-ii))), size(coeff_matrix, 1)) &
        !                   + 1, jj)
        !     end do 
        ! end do 

        do ii = 1, N 
            den(ii) = D**(N-ii)
        end do 
        do ii = 1, D**N 
            do jj = 1, N 
                psi(ii) = psi(ii) * coeff_matrix( mod((ii-1)/den(jj), D) + 1, jj)
            end do
        end do 

    end subroutine BuildSeparableState
            
    
    subroutine ReducedDensityMatrixB(rho, D, rho_B) 

        ! computes the reduced density matrix for N=2

        integer :: D, jj, kk, ll
        complex*16, dimension(:,:) :: rho 
        ! the reduced density matrix has dimension D^{N-1} x D^{N-1}
        complex*16, dimension(size(rho,1)/D, size(rho,2)/D) :: rho_B


        rho_B = 0.d0
        do jj = 1, size(rho,1)/D  
            do kk = jj, size(rho,1)/D 
                do ll = 1, size(rho,1)/D 
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
        ! the reduced density matrix has dimension D^{N-1} x D^{N-1}
        complex*16, dimension(size(rho,1)/D, size(rho,2)/D) :: rho_A

        rho_A = 0.d0
        do jj = 1, size(rho,1)/D  
            do kk = jj, size(rho,1)/D 
                do ll = 1, size(rho,1)/D
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

end module Utilities