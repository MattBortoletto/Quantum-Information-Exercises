module debug

    implicit none

    interface checkpoints
        module procedure check_int2, check_int4, check_real4, check_real8, check_cmplx8, check_cmplx16
        module procedure check_cmplxmat8, check_cmplxmat16
    end interface checkpoints

contains

    subroutine check_int2(debug, variable, to_print)
        logical :: debug 
        integer*2 :: variable
        character(*) :: to_print
        if (debug) then
            print *, "[Checkpoint]", to_print
        end if
    end subroutine check_int2

    subroutine check_int4(debug, variable, to_print)
        logical :: debug 
        integer*4 :: variable
        character(*) :: to_print
        if (debug) then
            print *, "[Checkpoint]", to_print
        end if
    end subroutine check_int4

    subroutine check_real4(debug, variable, to_print)
        logical :: debug 
        real*4 :: variable
        character(*) :: to_print
        if (debug) then
            print *, "[Checkpoint]", to_print
        end if
    end subroutine check_real4

    subroutine check_real8(debug, variable, to_print)
        logical :: debug 
        real*8 :: variable
        character(*) :: to_print
        if (debug) then
            print *, "[Checkpoint]", to_print
        end if
    end subroutine check_real8

    subroutine check_cmplx8(debug, variable, to_print)
        logical :: debug 
        complex*8 :: variable
        character(*) :: to_print
        if (debug) then
            print *, "[Checkpoint]", to_print
        end if
    end subroutine check_cmplx8

    subroutine check_cmplx16(debug, variable, to_print)
        logical :: debug 
        complex*16 :: variable
        character(*) :: to_print
        if (debug) then
            print *, "[Checkpoint]", to_print
        end if
    end subroutine check_cmplx16

    subroutine check_cmplxmat8(debug, variable, to_print)
        logical :: debug 
        complex*8, dimension(:, :) :: variable
        character(*) :: to_print
        if (debug) then
            print *, "[Checkpoint]", to_print
        end if
    end subroutine check_cmplxmat8

    subroutine check_cmplxmat16(debug, variable, to_print)
        logical :: debug 
        complex*16, dimension(:, :) :: variable 
        character(*) :: to_print
        if (debug) then
            print *, "[Checkpoint]", to_print
        end if
    end subroutine check_cmplxmat16
    