module Utilities

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

        write (str, "(F10.3)") k
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

end module Utilities