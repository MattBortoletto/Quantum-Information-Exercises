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
        call dfftw_destroy_plan 

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
        call dfftw_plan_dft_1d(n, size(array), array, ift, FFTW_FORWARD, FFTW_ESTIMATE)
        ! compute the actual transform 
        call dfftw_execute_dft(n, array, ift) 
        ! deallocate
        call dfftw_destroy_plan 

        ift = ift / sqrt(real(size(array)))

        return 

    end function InverseFourierTransform


    subroutine psiTimeEvol(psi_t, x_grid, p_grid, t_len, dt, omega, m)

        ! this subroutine computes the time evolution of the psi subject to 
        ! the "shifted" harmonic potential 

        complex*16, dimension(:,:) :: psi_t
        complex*16, dimension(:), allocatable :: evol_op_V, evol_op_T 
        !real*8, dimension(:,:), allocatable :: V_t 
        real*8, dimension(:) :: x_grid, p_grid 
        real*8 :: t_i, dt, omega, m 
        integer :: t_len, ii

        allocate(evol_op_V(t_len))
        allocate(evol_op_T(t_len))
        !allocate(V_t(N, t_len))

        evol_op_T = exp(dcmplx(0.d0, -dt*(p_grid)**2/(2*m)))
        
        do ii = 2, t_len

            t_i = ii * dt 
            evol_op_V = exp(dcmplx(0.d0,dt*m*(omega**2)*(x_grid-t_i)**2)/4)
            
            ! apply the position space operator 
            psi_t(:, ii) = psi_t(:, ii) * evol_op_V
            
            ! flip to momentum space using the FT 
            psi_t(:, ii) = FourierTransform(psi_t(:, ii))

            ! apply the momentum space operator
            psi_t(:, ii) = psi_t(:, ii) * evol_op_T

            ! flip to position space using the IFT
            psi_t(:, ii) = InverseFourierTransform(psi_t(:, ii))

            ! apply again the position space operator 
            psi_t(:, ii) = psi_t(:, ii) * evol_op_V

            ! checking if the norm is preserved (up to some finite precision)
            ! write(int_to_str,'(i3.3)') ii
            ! call debug(deb,psi_t(:,1),abs(sum(abs(psi_t(:,1)*sqrt(step))**2)&
            !     -1).ge.1d-4,&
            !     'no norm after '//trim(int_to_str)//' unitary(?) operator')
            
        end do

        return 

    end subroutine psiTimeEvol

end module TimeEvolution