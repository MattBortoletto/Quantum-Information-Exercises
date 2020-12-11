program qim

    use Utilities
    use debugger 
    use DensityMatrices
    use Ising  

    implicit none 

    integer :: N, kk, info, ii, ll  
    !real*8 :: lambda 
    complex*16, dimension(:,:), allocatable :: isingH 
    real*8, dimension(:), allocatable :: eig, lambda 
    !real*8 :: start, end 

    allocate(lambda(21))
    kk = 4
    lambda = (/((ii*0.15), ii=0,20)/) 

    do N = 2, 10

        open(unit=73, file='eig_'//trim(str_i(N))//'.txt', action="write")

        do ll = 1, size(lambda)

            print *, "---- N:", N, "--------------- lambda:", lambda(ll), "----" 

            !print *, "init"

            !call cpu_time(start)

            ! ---- initialize the Hamiltonian ----
            allocate(isingH(2**N, 2**N), stat=info)
            call checkpoint(debug=(info.ne.0), message="Allocation failed!", &
                            end_program=.true.)
            call IsingHamiltonian(isingH, N, lambda(ll))
            ! ------------------------------------

            !print *, "diag" 

            ! ---- diagonalize -------------------
            allocate(eig(2**N))
            call checkpoint(debug=(info.ne.0), message="Allocation failed!", &
                            end_program=.true.)
            call Diagonalize(isingH, eig, info)
            call checkpoint(debug=(info.ne.0), message="Diagonalization failed!", &
                            end_program=.true.)
            ! ------------------------------------

            !call cpu_time(end)
            !print *, end - start 

            !print *, "finish"

            ! ---- save the first kk eigenvalues divided by N (energy density) ---- 
            do ii = 1, size(eig(:kk)) 
                write(73, *) lambda(ll), eig(:kk) / (N-1)
            end do 
            ! ------------------------------------

            ! ---- deallocate memory -------------
            deallocate(eig, isingH)
            ! ------------------------------------

        end do 

        close(73)

    end do 

end program qim
