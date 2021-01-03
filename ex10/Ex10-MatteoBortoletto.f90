program rsrg_ising

    use Utilities
    use debugger 
    use DensityMatrices
    use Ising  

    implicit none 

    ! N: number of subsystems
    ! info: stat flag
    ! ii, ll: variables to loop
    ! niter: number of iterations for the RSRG
    integer :: N, info, ii, ll, niter
    ! Ising Hamiltonian
    real*8, dimension(:,:), allocatable :: isingH 
    ! vector of field strenght values
    real*8, dimension(:), allocatable :: lambda 
    ! ground state value
    real*8 :: gs 

    allocate(lambda(21))

    niter = 100
    lambda = (/((ii*0.15), ii=0,20)/) 

    do N = 2, 5

        open(unit=73, file='gs_N'//trim(str_i(N))//'.txt', action="write")

        do ll = 1, size(lambda)

            print *, "---- N:", N, "--------------- lambda:", lambda(ll) 

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
            write(73, *) lambda(ll), gs / (N*2.0**(niter+1))
            ! ---------------------------------------------------------------------

            ! ---- deallocate memory ----------------------------------------------
            deallocate(isingH)
            ! ---------------------------------------------------------------------

        end do 

        close(73)

    end do 

end program rsrg_ising 
