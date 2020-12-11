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
            write(73, *) vect(ii)
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

end module Utilities