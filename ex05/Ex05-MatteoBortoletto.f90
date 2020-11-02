program eigenproblem 

    ! EX01 - EIGENPROBLEM
    ! main logic:
    ! - ask the user to insert the dimension n of the matrix
    ! - fill it with random numbers
    ! - call the subroutine to compute the eigenvalues 
    ! - compute the s_i (create a subroutine/function)
    ! - optional point

    implicit none

    ! dimension of the matrix
    integer n 
    ! matrix
    real*8, dimension(:,:), allocatable :: m
    ! array to store the eigenvalues 
    real*8, dimension(:), allocatable :: eig 

    ! ask the user to enter the dimension of the matrix
    print *, "Please enter the dimension of the matrix: "
    read *, n 

    ! allocate the memory
    allocate(m(n, n))

    ! call the subroutine to compute the eigenvalues 
    
    
    ! save the results in a text file 
    ! open(10, file='results.txt', status='replace')
    ! write(10, *) 'Eigenvalues:'
    ! do ii = 1, n
    !     write(10, 10) ii, eig(ii)
    !     10 format(I3, '   ', f14.8)
    ! end do
    ! write(10, *)
    ! write(10, *) 'Eigenvectors:'
    ! do ii = 1, n
    !     write(10, 20) ii, m(:, ii)
    !     20 format(i3, '   ', 10f14.8)
    ! end do
    ! write(10, *)
    
    ! close(10)
    
    deallocate(m)
    deallocate(eig)

    stop 
    
end program eigenproblem
