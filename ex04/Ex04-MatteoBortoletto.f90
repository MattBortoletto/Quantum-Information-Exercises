module debugger 

    ! module for debugging
    ! It contains a subroutine that can be used as checkpoint for debugging. 

    implicit none

contains

    subroutine checkpoint(debug, variable, message, end_program)

        implicit none

        ! debug: logical variable that enables debugging if it is .true.
        ! end_program: logical variable that stops the program if is .true.
        logical, intent(in) :: debug, end_program
        ! optional variable to be printed 
        class(*), intent(in), optional :: variable
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
            print *, "----------------------------------------------------------------------"
            ! if end_program = true then in case of error the program will stop
            if (end_program) then 
                stop
            end if 
        end if 

    end subroutine checkpoint

end module debugger 


module mult 

    ! module for matrix multiplication 

    implicit none 

    ! matrices (m1 x m2 = m3)
    real*8, dimension(:,:), allocatable :: m1, m2, m3
    ! indices for loop and dimensions of the matrices
    integer :: ii, jj, kk, ll, mm, nn, pp, rr, nrows1, ncols1, nrows2, ncols2
    ! timing variables
    real*8 :: time_start1, time_stop1, time_start2, time_stop2, time_start3, time_stop3, time1, time2, time3
    ! variable to store the user's aswer (y/n) when asked if he/she wants to save the results of the multiplication
    character(1) :: save_results

contains

    function mult1(m_1, m_2) result(m_3) 
        ! non optimized function for matrix multiplication

        ! pre: real*8 matrix multiplication.
        !      the number of columns of m1 must be equal to the number of rows of m2.
        ! post: a real*8 matrix m3 which elements are the dot product of the i-th row 
        !       of m1 and the j-th column of m2 is returned.

        implicit none 

        real*8, dimension(nrows1, ncols1) :: m_1 
        real*8, dimension(ncols1, ncols2) :: m_2
        real*8, dimension(nrows1, ncols2) :: m_3

        m_3 = 0.
        do ii = 1, nrows1
            do jj = 1, ncols2 
                do kk = 1, ncols1
                    m_3(ii, jj) = m_3(ii, jj) + m_1(ii, kk) * m_2(kk, jj)
                end do
            end do
        end do

    end function mult1

    function mult2(m_1, m_2) result(m_3) 
        ! optimized function for matrix multiplication

        ! pre: real*8 matrix multiplication.
        !      The number of columns of m1 must be equal to the number of rows of m2.
        ! post: a real*8 matrix m3 which elements are the dot product of the i-th row 
        !       of m1 and the j-th column of m2 is returned.

        implicit none

        real*8, dimension(nrows1, ncols1) :: m_1 
        real*8, dimension(ncols1, ncols2) :: m_2
        real*8, dimension(nrows1, ncols2) :: m_3

        m_3 = 0.
        ! Since in Fortran matrices are stored in memory as columns, 
        ! we can set the inner loop to be over consecutive elements 
        ! in memory in this way:
        do jj = 1, ncols2
            do kk = 1, ncols1
                do ii = 1, nrows1
                    m_3(ii, jj) = m_3(ii, jj) + m_1(ii, kk) * m_2(kk, jj)
                end do
            end do
        end do

    end function mult2

    subroutine LoadDimensions(filename, dimensions)
        ! this subroutine is used to load the dimension of the two 
        ! matrices from a txt files 

        implicit none

        ! name of the file
        character(*), intent(in) :: filename
        ! array to store the dimensions
        integer, dimension(4) :: dimensions
        ! value for the iostat specifier
        integer :: stat 

        open(unit = 3, file = filename, action = "read", status = "old")
        stat = 0
        ii = 1
        do while(stat == 0)
            read(3, *, iostat = stat) dimensions(ii)
            ii = ii + 1
        end do

        close(3)

    end subroutine LoadDimensions

end module mult



program MyMatrixMultiplication

    ! module for matrix multiplication
    use mult
    ! module for debugging 
    use debugger

    implicit none

    integer, dimension(4) :: dims

    ! load the dimension of the matrices
    call LoadDimensions('matrix_dimensions.txt', dims)
    
    ! assign the dimensions
    nrows1 = dims(1)
    ncols1 = dims(2)
    nrows2 = dims(3)
    ncols2 = dims(4)

    ! ask to enter the dimension until nrows and ncols are greater or equal than 1 and nrows2=ncols1
    do while ((nrows2 .ne. ncols1) .or. (((nrows1 .lt. 1) .or. (ncols1 .lt. 1)) .or. ((nrows2 .lt. 1) .or. (ncols2 .lt. 1))))
        ! display two different messages according to the condition which is not satisfied

        ! if nrows2=ncols1 warn the user that the number of rows of the second matrix must 
        ! be equal to the number of columns of the first matrix
        if (nrows2 .ne. ncols1) then 
            call checkpoint(debug = .true., message = "The number of columns in the first matrix is not equal &
                                                        &to the number of rows in the second matrix.", end_program = .false.)
            print *, "Please modify the text file."
            stop
        end if 

        ! if some of the dimensions is less than 1 warn the user and ask to enter again
        if (((nrows1 .lt. 1) .or. (ncols1 .lt. 1)) .or. ((nrows2 .lt. 1) .or. (ncols2 .lt. 1))) then 
            call checkpoint(debug = .true., message = "Non valid dimension! The number of rows and columns must &
                                                        &be greater or equal than 1.", end_program = .false.)
            print *, "Please modify the text file."
            stop 
        end if 
    end do

    ! allocate the memory
    allocate(m1(nrows1, ncols1))
    allocate(m2(nrows2, ncols2))
    allocate(m3(nrows1, ncols2))

    ! fill the matrices with random integer numbers
    call random_number(m1)
    call random_number(m2)
    m1 = int(m1*10)
    m2 = int(m2*10)

    ! open a file to save the results of the computations
    ! warning: if the matrices are big, the txt file can be very large (~100 MB)
    print *, 'Do you want to save the results in a txt file (in case of big matrices it can be very large, ~100 MB) [y/n]?'
    read *, save_results
    if (save_results == 'y') then 
        open(unit = 2, file = 'results.txt', action = "write", status = "replace") 
        write(2, *) "m1:"
        do ll = 1, nrows1
            write(2, *) real(m1(ll, :)) 
        end do
        write(2, *) " "
        write(2, *) "m2:"
        do mm = 1, ncols1
            write(2, *) real(m2(mm, :))
        end do
    end if

    ! non optimized function
    call cpu_time(time_start1)
    m3 = mult1(m1, m2)
    call cpu_time(time_stop1)
    time1 = time_stop1 - time_start1
    print *, "Standard loop time = ", time1

    ! write the results in the txt file
    if (save_results == 'y') then
        write(2, *) " "
        write(2, *) "m3 with mult1:"
        do nn = 1, nrows1
            write(2, *) real(m3(nn, :))
        end do
    end if

    ! optimized function
    call cpu_time(time_start2)
    m3 = mult2(m1, m2)
    call cpu_time(time_stop2)
    time2 = time_stop2 - time_start2
    print *, "Optimized loop time = ", time2

    ! check if the optimized function is faster
    if (time2 < time1) then 
        print *, "The optimized loop is faster!"
    else 
        print *, "The optimized loop is not faster!"
    end if 

    ! write the results in the txt file
    if (save_results == 'y') then
        write(2, *) " "
        write(2, *) "m3 with mult2:"
        do pp = 1, nrows1
            write(2, *) real(m3(pp, :))
        end do
    end if

    ! instrinsic function
    call cpu_time(time_start3)
    m3 = matmul(m1, m2)
    call cpu_time(time_stop3)
    time3 = time_stop3 - time_start3
    print *, "Matmul time = ", time3

    ! write the results in the txt file
    if (save_results == 'y') then
        write(2, *) " "
        write(2, *) "m3 with matmul:"
        do rr = 1, nrows1
            write(2, *) real(m3(rr, :))
        end do
    end if 

    ! close the txt file
    close(2)

    deallocate(m1, m2, m3)

    stop 

end program MyMatrixMultiplication