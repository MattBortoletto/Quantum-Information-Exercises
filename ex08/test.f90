program prova

    implicit none 

    integer :: D, N, ii, jj 
    integer, dimension(:), allocatable :: den 

    D = 2
    N = 5

    allocate(den(N))

    do ii = 1, N 
        den(ii) = D**(N-ii)
        print *, "den", ii, den(ii)
    end do

    do ii = 1, D**N 
        do jj = 1, N 
            print *, "mod", mod((ii-1)/den(jj), D) 
        end do
    end do 

end program