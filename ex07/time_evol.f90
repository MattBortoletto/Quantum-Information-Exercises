module TimeEvolution 

    use, intrinsic :: iso_c_binding

    implicit none 

    include 'fftw3.f03'

    contains 

    function FourierTransform(array) result(ft) 

        ! this function computed the Fourier transfor 
        ! using the FFTW package

        complex*16, dimension(:) :: array 
        complex*16, dimension(size(array)) :: ft 
        integer*8 :: n

        ! create a plan 
        call dfftw_plan_dft_1d(n, size(array), array, ft, FFTW_FORWARD, FFTW_ESTIMATE)
        ! compute the actual transform 
        call dfftw_execute_dft(n, array, ft) 
        ! deallocate
        call dfftw_destroy_plan(n)

        ft = ft / sqrt(real(size(array)))

        return 

    end function FourierTransform 


    function InverseFourierTransform(array) result(ift)

        ! this function computed the Fourier transfor 
        ! using the FFTW package

        complex*16, dimension(:) :: array 
        complex*16, dimension(size(array)) :: ift 
        integer*8 :: n

        ! create a plan 
        call dfftw_plan_dft_1d(n, size(array), array, ift, FFTW_BACKWARD, FFTW_ESTIMATE)

        ! compute the actual transform 
        call dfftw_execute_dft(n, array, ift) 

        ! deallocate
        call dfftw_destroy_plan(n) 

        ift = ift / sqrt(real(size(array))) 

        return 

    end function InverseFourierTransform


    subroutine psiTimeEvol(psi_t, prob_t, x_grid, p_grid, t_len, dt, omega, m, hbar, T)

        ! this subroutine computes the time evolution of the psi subject to 
        ! the "shifted" harmonic potential 

        complex*16, dimension(:,:) :: psi_t
        real*8, dimension(:,:) :: prob_t 
        complex*16, dimension(:), allocatable :: tmp 
        complex*16, dimension(:), allocatable :: evol_op_V, evol_op_T 
        real*8, dimension(:) :: x_grid, p_grid 
        real*8 :: t_i, dt, omega, m, hbar, T 
        integer :: t_len, ii

        allocate(evol_op_V(t_len))
        allocate(evol_op_T(t_len))
        allocate(tmp(size(psi_t, 1)))

        tmp = psi_t(:, 1)

        ! evol_op_T = exp(dcmplx(0.d0, -dt*(p_grid)**2/(2*m*hbar)))
        evol_op_T = exp( -dt * dcmplx(0.d0,(p_grid)**2) / (2*m*hbar) )

        do ii = 1, t_len-1

            t_i = ii * dt 
            ! evol_op_V = exp(-dt*m * dcmplx(0.d0,(omega**2)*(x_grid-t_i/T)**2) / (4*hbar) ) 
            evol_op_V = exp(-dt*m * dcmplx(0.d0,(omega**2)*(x_grid-t_i)**2) / (4*hbar) ) 

            ! apply the position space operator 
            tmp = tmp * evol_op_V

            ! flip to momentum space using the FT 
            tmp = FourierTransform(tmp)

            ! apply the momentum space operator
            tmp = tmp * evol_op_T

            ! flip to position space using the IFT
            tmp = InverseFourierTransform(tmp)

            ! apply again the position space operator 
            tmp = tmp * evol_op_V

            psi_t(:, ii+1) = tmp 
            
            prob_t(:, ii+1) = abs(tmp)**2
            
        end do

        deallocate(tmp) 

        return 

    end subroutine psiTimeEvol


    subroutine update(psi, N, dx, q0, xmax, dt, dp)

        integer :: ii, N 
        complex*16, dimension(:) :: psi
        complex*16, dimension(:), allocatable :: psi_fft 
        real*8 :: dx, q0, xmax, dt, dp 
        integer*8 :: plan 

        allocate(psi_fft(size(psi)))

        do ii = 1, N 
            psi(ii) = psi(ii)*exp(-0.5*cmplx(0, (ii*dx-q0-xmax)**2)*dt)
        end do 

        call dfftw_plan_dft_1d(plan, N, psi, psi_fft, FFTW_FORWARD, FFTW_ESTIMATE)
        call dfftw_execute_dft(plan, psi, psi_fft) 
        call dfftw_destroy_plan(plan)

        do ii = 1, N/2
            psi_fft(ii) = psi_fft(ii)*exp(-0.5*cmplx(0, ((ii-1)*dp)**2)*dt)
        end do 

        do ii = N/2+1, N 
            psi_fft(ii) = psi_fft(ii)*exp(-0.5*cmplx(0, ((ii-1-N)*dp)**2)*dt)
        end do 

        call dfftw_plan_dft_1d(plan, N, psi_fft, psi, FFTW_BACKWARD, FFTW_ESTIMATE)
        call dfftw_execute_dft(plan, psi_fft, psi) 
        call dfftw_destroy_plan(plan)

        do ii = 1, N 
            psi(ii) = psi(ii)*exp(-0.5*cmplx(0, (ii*dx-q0-xmax)**2)*dt)/N 
        end do
        
    end subroutine update 


end module TimeEvolution