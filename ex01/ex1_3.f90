module mult 

    real*8, dimension(:,:), allocatable :: m1, m2, m3
    integer*2 :: ii, jj, kk, nrows1, ncols1, ncols2, ll, mm, nn, pp, rr
    real*8 :: time_start1, time_stop1, time_start2, time_stop2, time_start3, time_stop3, time1, time2, time3
    character (len = 1) :: save_results

contains

function mult1(m_1, m_2) result(m_3) ! non optimized function

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

function mult2(m_1, m_2) result(m_3) ! optimized function

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

end module mult


program MyMatrixMultiplication

    use mult
    implicit none
    
    nrows1 = 1200
    ncols1 = 1300
    ncols2 = 1100

    allocate(m1(nrows1, ncols1))
    allocate(m2(ncols1, ncols2))
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