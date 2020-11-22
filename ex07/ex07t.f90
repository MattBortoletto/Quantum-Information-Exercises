module DiagNorm

    use, intrinsic :: iso_c_binding

    implicit none

    include 'fftw3.f03'

    contains

    subroutine Diagonalize(d, ud, ev)

        real*8 :: abstol
        integer :: m, info 
        real*8, dimension(:), allocatable :: d, ud, eig, work
        integer, dimension(:), allocatable :: iwork, ifail 
        real*8, dimension(size(d), 1) :: ev 

        allocate(eig(size(d)))
        allocate(work(5*size(d)))
        allocate(iwork(5*size(d)))
        allocate(ifail(size(d)))

        abstol = -1

        call dstevx('v', 'i', size(d), d, ud, 0d0, 10d0, 1, 1, abstol, &
                    m, eig, ev, size(d) ,work, iwork, ifail, info)
        
        if (info /= 0) then
            print *,  "diagonalization failed!"
        end if 

        if (abs(eig(1)-0.5)>0.01) then 
            print *, "wrong ground state energy!"
        end if 

        deallocate(d, ud, eig, work, iwork, ifail)

    end subroutine Diagonalize

    subroutine NormalizePsi(psi, dx) 

        complex*16, dimension(:) :: psi 
        real*8 :: dx, norm 
        integer :: ii 

        norm = 0.0

        do ii = 1, size(psi) 
            norm = norm + dx*(real(psi(ii))**2 + aimag(psi(ii))**2)
        end do 

        psi = psi / sqrt(norm) 

        do ii = 1, size(psi)
            norm = norm + dx*(real(psi(ii))**2 + aimag(psi(ii))**2)
        end do

        norm = 0.0

        do ii = 1, size(psi) 
            norm = norm + dx*(real(psi(ii))**2 + aimag(psi(ii))**2)
        end do 

        if (abs(norm-1.0) .gt. 0.00001) then 
            print *, "Normalization error!"
            stop 
        end if

    end subroutine NormalizePsi

    subroutine TimeEvol(psi_x, psi_p, T, dt, dx, dp, L)

        integer :: N_t, tt, ii
        real*8 :: T, dt, min, dx, dp, L, mean 
        complex*16 :: imm 
        integer :: plan
        complex*16, dimension(:) :: psi_x, psi_p
        

        N_t = int(T/dt)

        imm = cmplx(0.0, 1.0)

        ! do tt = 1, N_t

        !     print *, tt 

        !     min = L/2 + tt*dt/T 

        !     do ii = 1, size(psi_x) 
        !         psi_x(ii) = psi_x(ii) * exp(-imm*dt*0.25*(((ii-1)*dx - min)**2))
        !     end do 

        !     call dfftw_plan_dft_1d(plan, size(psi_x), psi_x, psi_p, &
        !                         FFTW_FORWARD, FFTW_ESTIMATE)
        !     call dfftw_execute_dft(plan, psi_x, psi_p)
        !     call dfftw_destroy_plan(plan)

        !     call Normalize(psi_p, dp)

        !     do ii = 2, (size(psi_x)+1)/2
        !         psi_p(ii) = psi_p(ii) * exp(-imm*dt*0.5*(((ii-1)*dp)**2))
        !         psi_p(size(psi_x)+2-ii) = psi_p(size(psi_x)+2-ii) * exp(-imm*dt*0.5*(((ii-1)*dp)**2))
        !     end do 

        !     psi_p((size(psi_x)+2)/2) = psi_p((size(psi_x)+2)/2) * exp(-imm*dt*0.5*(((atan(1.0)/dx)**2)))

        !     call dfftw_plan_dft_1d(plan, size(psi_p), psi_p, psi_x, &
        !                         FFTW_BACKWARD, FFTW_ESTIMATE)
        !     call dfftw_execute_dft(plan, psi_p, psi_x)
        !     call dfftw_destroy_plan(plan)

        !     call Normalize(psi_x, dx)

        !     do ii = 1, size(psi_x)
        !         psi_x(ii) = psi_x(ii) * exp(-imm*dt*0.25*(((ii-1)*dx - min)**2))
        !     end do

        !     call Normalize(psi_x, dx)

        !     mean = 0.0
        !     do ii = 1, size(psi_x)
        !         mean = mean + (dx**2)*((real(psi_x(ii)))**2+(aimag(psi_x(ii)))**2)*(ii-1)
        !     end do 

        !     open(73, file="mean_pos.txt", access="append")
        !     write(73,*) min, mean
        !     close(73)

        ! end do 

        do tt = 1, N_t

            print *, tt 

            min = L/2 + tt*dt/T

            do ii = 1, size(psi_x)
                psi_x(ii)=exp(-imm*dt*0.25*(((ii-1)*dx-min)**(2)))*psi_x(ii)
            end do

            call dfftw_plan_dft_1d(plan, size(psi_x), psi_x, psi_p, FFTW_FORWARD, FFTW_ESTIMATE)
            call dfftw_execute_dft(plan, psi_x, psi_p)
            call dfftw_destroy_plan(plan)

            call NormalizePsi(psi_p, dp)

            do ii = 2, (size(psi_p)+1)/2
                psi_p(ii)=exp(-imm*dt*0.5*(((ii-1)*dp)**(2)))*psi_p(ii)
                psi_p(size(psi_p)+2-ii)=exp(-imm*dt*0.5*(((ii-1)*dp)**(2)))*psi_p(size(psi_p)+2-ii)
            end do
            psi_p((size(psi_p)+2)/2)=exp(-imm*dt*0.5*((atan(1.0)/dx)**(2)))*psi_p((size(psi_p)+2)/2)

            call dfftw_plan_dft_1d(plan, size(psi_p), psi_p, psi_x, FFTW_BACKWARD, FFTW_ESTIMATE)
            call dfftw_execute_dft(plan, psi_p, psi_x)
            call dfftw_destroy_plan(plan)

            call NormalizePsi(psi_x, dx)

            do ii = 1, size(psi_x)
                psi_x(ii)=exp(-imm*0.25*dt*(((ii-1)*dx-min)**(2)))*psi_x(ii)
            end do

            call NormalizePsi(psi_x, dx)

            mean=0.0

            do ii = 1, size(psi_x)
                mean = mean+dx*dx*((real(psi_x(ii)))**(2)+(aimag(psi_x(ii)))**(2))*(ii-1)
            end do

            open(13, file="mean_pos.txt", ACCESS="append")
            write(13,*) min, mean
            close(13)

        end do

    end subroutine TimeEvol

end module DiagNorm


program time_evol

    use DiagNorm

    implicit none 

    real*8 :: dx, dp, dt, T, L, min
    real*8, dimension(:,:), allocatable :: ev 
    integer :: bins, ii
    complex*16, dimension(:), allocatable :: psi_x, psi_p
    real*8, dimension(:), allocatable :: d, ud
    
    dx = 0.001
    dt = 0.1
    T = 1000
    L = 50

    bins = int(L/dx) + 1
    min = L/2
    dp = atan(1.0)/L 

    allocate(d(bins))
    allocate(ud(bins-1))
    allocate(ev(bins, 1))

    do ii = 1, size(d) 
        d(ii) = 2.0/(4.0*dx*dx) + ((ii-1)*dx - min)**2
        if (ii == size(d)) then
            exit 
        end if 
        ud(ii) = - 1.0 / (4.0*dx*dx)
    end do 

    call Diagonalize(d, ud, ev)

    allocate(psi_x(bins), psi_p(bins))

    do ii = 1, bins 
        psi_x(ii) = cmplx(ev(ii, 1), 0.0)
    end do 

    call Normalize(psi_x, dx) 

    call TimeEvol(psi_x, psi_p, T, dt, dx, dp, L)

    open(42, file="prob_ev.txt")
    do ii = 1, bins
        write(42, *) (ii-1)*dx, real(psi_x(ii))**2 + aimag(psi_x(ii))**2
    end do 
    close(42)

end program time_evol