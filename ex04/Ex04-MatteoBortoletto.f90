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
    !character(1) :: save_results

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
        integer :: iostat 

        open(unit = 10, file = filename, status = "old")
        ii = 1
        do while(ii < 5)
            read(10, *, iostat=iostat) dimensions(ii)
            if(iostat < 0) then
                write(6,'(A)') "Warning: File containts less than 4 entries."
                exit
            else if(iostat > 0) then
                write(6,'(A)') 'Error reading file.'
                stop
            end if
            ii = ii + 1
        end do

        close(10)

        ! -------------------------------------------------------------------------------

        ! open(unit = 3, file = filename, action = "read", status = "old")
        ! ! if (stat .ne. 0) then 
        ! !     print *, "Cannot open file"
        ! !     stop
        ! ! end if 
        ! do ii = 1, 4
        !     ! read(3, *, iostat = stat) dimensions(ii)
        !     read(3, *) dimensions(ii)
        ! end do

        ! close(3)

    end subroutine LoadDimensions

    subroutine StoreResult(filename, var)
        ! this subroutine is used to store some results in a txt file
        ! in particular, it will be used to store the CPU times for the
        ! different matrix multiplication methods

        implicit none

        !name of the file
        character(*), intent(in) :: filename 
        ! variable to store the time 
        real*8, intent(in) :: var 
        ! logical variable to see if the file 'filename' altrady exists
        logical :: exists 
        
        ! check if the file exists
        inquire(file=filename, exist=exists)
        ! if the file exists then open it, otherwise create it 
        if (exists) then 
            open(unit = 3, file = filename, action = "write", status = "old", position = "append")
        else 
            open(unit = 3, file = filename, action = "write", status = "new")
        end if 
        
        ! write the result
        write(3, *) nrows1, var
        ! close the file 
        close(3)

    end subroutine StoreResult

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

    ! ----------------------------
    print *, 'fortran dims:', dims
    ! ----------------------------

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

    ! non optimized function
    call cpu_time(time_start1)
    m3 = mult1(m1, m2)
    call cpu_time(time_stop1)
    time1 = time_stop1 - time_start1
    !print *, "Standard loop time = ", time1
    call StoreResult("not-optimized.txt", time1)

    ! optimized function
    call cpu_time(time_start2)
    m3 = mult2(m1, m2)
    call cpu_time(time_stop2)
    time2 = time_stop2 - time_start2
    !print *, "Optimized loop time = ", time2
    call StoreResult("optimized.txt", time2)

    ! instrinsic function
    call cpu_time(time_start3)
    m3 = matmul(m1, m2)
    call cpu_time(time_stop3)
    time3 = time_stop3 - time_start3
    !print *, "Matmul time = ", time3
    call StoreResult("matmul.txt", time3)

    deallocate(m1, m2, m3)

    stop 

end program MyMatrixMultiplication