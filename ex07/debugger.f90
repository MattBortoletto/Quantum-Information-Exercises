module debugger 

    ! module for debugging
    ! It contains a subroutine that can be used as checkpoint for debugging. 

    implicit none

contains

    subroutine checkpoint(debug, variable, array_variable, matrix_variable, message, end_program)

        implicit none

        ! debug: logical variable that enables debugging if it is .true.
        ! end_program: logical variable that stops the program if is .true.
        logical, intent(in) :: debug, end_program
        ! optional variable to be printed 
        class(*), intent(in), optional :: variable
        ! optional array to be printed
        class(*), dimension(:), intent(in), optional :: array_variable
        ! optional matrix to be printed
        class(*), dimension(:,:), intent(in), optional :: matrix_variable
        ! optional message to be printed 
        character(*), optional :: message 

        ! if debug = true then the message will be printed and in case of presence of a variable
        ! it will also print the it
        if (debug) then
            print *, "[Debugger] -----------------------------------------------------------"
            print *, message
            if (present(variable)) then
                select type(variable) 
                    type is (integer(2))
                        print *, variable
                    type is (integer(4))
                        print *, variable
                    type is (real(4))
                        print *, variable
                    type is (real(8))
                        print *, variable
                    type is (complex(8))
                        print *, variable
                    type is (complex(16))
                        print *, variable
                    type is (logical)
                        print *, variable
                end select
            end if 
            if (present(array_variable)) then
                select type(array_variable) 
                    type is (integer(2))
                        print *, array_variable
                    type is (integer(4))
                        print *, array_variable
                    type is (real(4))
                        print *, array_variable
                    type is (real(8))
                        print *, array_variable
                    type is (complex(8))
                        print *, array_variable
                    type is (complex(16))
                        print *, array_variable
                end select 
            end if 
            if (present(matrix_variable)) then
                select type(matrix_variable) 
                    type is (integer(2))
                        print *, matrix_variable
                    type is (integer(4))
                        print *, matrix_variable
                    type is (real(4))
                        print *, matrix_variable
                    type is (real(8))
                        print *, matrix_variable
                    type is (complex(8))
                        print *, matrix_variable
                    type is (complex(16))
                        print *, matrix_variable
                end select 
            end if 
            print *, "----------------------------------------------------------------------"
            ! if end_program = true then in case of error the program will stop
            if (end_program) then 
                stop
            end if 
        end if 

    end subroutine checkpoint

end module debugger 