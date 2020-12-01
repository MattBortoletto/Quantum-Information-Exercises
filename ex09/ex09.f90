program qim

    use Utilities
    use debugger 
    use DensityMatrices
    use Ising  

    implicit none 

    integer :: N, kk, info 
    real*8 :: lambda 
    complex*16, dimension(:,:), allocatable :: isingH 
    real*8, dimension(:), allocatable :: eig
    real*8 :: start, end 

    !N = 2
    lambda = 1
    kk = 3

    do N = 2, 15

        print *, "-----", N, "------------------" 

        print *, "init"

        call cpu_time(start)

        ! ---- initialize the Hamiltonian ----
        allocate(isingH(2**N, 2**N), stat=info)
        call checkpoint(debug=(info.ne.0), message="Allocation failed!", &
                        end_program=.true.)
        call IsingHamiltonian(isingH, N, lambda)
        ! ------------------------------------

        print *, "diag" 

        ! ---- diagonalize -------------------
        allocate(eig(2**N))
        call checkpoint(debug=(info.ne.0), message="Allocation failed!", &
                        end_program=.true.)
        call Diagonalize(isingH, eig, info)
        call checkpoint(debug=(info.ne.0), message="Diagonalization failed!", &
                        end_program=.true.)
        ! ------------------------------------

        call cpu_time(end)

        print *, end - start 

        print *, "finish"

        ! ---- deallocate memory -------------
        deallocate(eig, isingH)
        ! ------------------------------------

    end do 

end program qim 
