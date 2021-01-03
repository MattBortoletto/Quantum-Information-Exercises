program rsrg_ising

    use Utilities
    use debugger 
    use DensityMatrices
    use Ising  

    implicit none 

    integer :: N, kk, info, ii, ll, niter
    real*8, dimension(:,:), allocatable :: isingH 
    real*8, dimension(:), allocatable :: lambda 
    real*8 :: gs 

    allocate(lambda(21))
    kk = 4
    niter = 100
    lambda = (/((ii*0.15), ii=0,20)/) 

    do N = 2, 4

        open(unit=73, file='gs_N'//trim(str_i(N))//'.txt', action="write")

        do ll = 1, size(lambda)

            print *, "---- N:", N, "--------------- lambda:", lambda(ll), "----" 

            ! ---- initialize the Hamiltonian -------------------------------------
            allocate(isingH(2**N, 2**N), stat=info)
            call checkpoint(debug=(info.ne.0), message="Allocation failed!", &
                            end_program=.true.)
            call IsingHamiltonian(isingH, N, lambda(ll))
            ! ---------------------------------------------------------------------

            ! ---- use Real Space Renormalization Group ---------------------------
            gs = RSRG(isingH, N, niter) 
            ! ---------------------------------------------------------------------

            ! ---- save the results -----------------------------------------------
            write(73, *) lambda(ll), gs / (n*2**(niter+1))
            ! ---------------------------------------------------------------------

            ! ---- deallocate memory ----------------------------------------------
            deallocate(isingH)
            ! ---------------------------------------------------------------------

        end do 

        close(73)

    end do 

end program rsrg_ising 
