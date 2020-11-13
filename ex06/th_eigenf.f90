! from https://sukhbinder.wordpress.com/hermite-polynomials-fortran-module/

MODULE hermite

CONTAINS

    RECURSIVE FUNCTION HermitePoly(n) RESULT(hp)

        REAL*8 hp(n + 1)
        REAL*8 hp1(n + 1), hp2(n + 1)

        IF (n .EQ. 0) THEN
            hp(1) = 1.0
            RETURN
        END IF

        IF (n .EQ. 1) THEN
            hp(1) = 2.0
            hp(2) = 0.0
        ELSE
            hp1(1:n + 1) = 0.0
            hp1(1:n) = 2.0*HermitePoly(n - 1)
            hp2(1:n + 1) = 0.0
            hp2(3:) = 2.0*(n - 1)*HermitePoly(n - 2)
            hp = hp1 - hp2
        END IF

    END FUNCTION

    FUNCTION evalHermitePoly(ix, n) RESULT(y)

        INTEGER n, ip
        REAL*8 ix(:), y(size(ix)), h(n + 1)

        k = size(ix)
        h = HermitePoly(n)

        y(1:k) = h(n + 1)
        ip = 1
        DO i = n, 1, -1
            DO j = 1, k
                y(j) = y(j) + h(i)*ix(j)**ip
            END DO
            ip = ip + 1
        END DO

    END FUNCTION

END MODULE HERMITE

! ----------------------------------------------------------------------------------

module Utilities

    implicit none

contains

    function str_i(k) result(str)

        ! converts an integer into string
        character(len=20) :: str
        integer, intent(in) :: k

        write (str, *) k
        str = adjustl(str)

        return

    end function str_i

    function str_r_d(k) result(str)

        ! converts a real into string (decimal notation)
        character(len=20) :: str
        real*8, intent(in) :: k

        write (str, "(F10.3)") k
        str = adjustl(str)

        return

    end function str_r_d

    function str_r_e(k) result(str)

        ! converts a real into string (exponential notation)
        character(len=20) :: str
        real*8, intent(in) :: k

        write (str, "(E10.1)") k
        str = adjustl(str)

        return

    end function str_r_e

    real function fact(n)

        integer, intent(in) :: n
        integer :: i

        if (n < 0) error stop 'factorial is singular for negative integers'

        fact = 1.0
        do i = 2, n
            fact = fact*i
        end do

    end function fact

    subroutine WriteEigenvalues(eig, filename)

        character(*) :: filename
        real*8, dimension(:) :: eig
        integer :: ii

        open (unit=73, file=filename, action="write", status="replace")
        do ii = 1, size(eig)
            write (73, *) eig(ii)
        end do
        close (73)

        return

    end subroutine WriteEigenvalues

    subroutine WriteEigenvectors(matr, grid_points, filename)

        character(*) :: filename
        real*8, dimension(:, :) :: matr
        real*8, dimension(:) :: grid_points
        integer :: ii

        open (unit=73, file=filename, action="write", status="replace")
        do ii = 1, size(matr, 1)
            write (73, *) grid_points(ii), matr(ii, :)
        end do
        close (73)

        return

    end subroutine WriteEigenvectors

    subroutine ComputeProb(H, prob)

        complex*16, dimension(:, :) :: H
        real*8, dimension(size(H, 1), size(H, 2)) :: prob
        integer :: ii 

        do ii = 1, size(H, 2)
            prob(:, ii) = abs(H(:, ii))**2
        end do

        return

    end subroutine ComputeProb

end module Utilities

! ----------------------------------------------------------------------------------

program hermite_polynomals

    use hermite
    use Utilities

    implicit none

    integer :: N, L, ii, jj, k
    real*8 :: dx, omega, m, hbar, pi
    complex*16, dimension(:, :), allocatable :: eigenf
    real*8, dimension(:), allocatable :: grid_points
    real*8, dimension(:, :), allocatable :: prob
    character(:), allocatable :: eigenf_filename, prob_filename
    character(1) :: which_param

    pi = 4.d0*datan(1.d0)

    print *, "Do you want to use use custom parameters of the default one (L=500, dx=0.01, omega=5, m=hbar=1)? [c/d]"
    read *, which_param
    if (which_param == "d") then 
        ! ---- default values ----
        L = 500
        dx = 0.01
        omega = 5
        m = 1.0
        hbar = 1.0
        ! ------------------------
    else if (which_param == "c") then 
        print *, "Please enter L, dx, omega, m and hbar:"
        read *, L, dx, omega, m, hbar 
    else 
        print *, "Invalid input."
        stop 
    end if 

    N = L*2 + 1

    print *, "Please enter the number of eigenfunctions you want to generate:"
    read *, k

    allocate (grid_points(N))
    do jj = 1, N
        ! x_min = -L*dx
        grid_points(jj) = -L*dx + dx*(jj - 1)
    end do

    allocate (eigenf(N, k))

    do ii = 1, k
        eigenf(:, ii) = (1/(sqrt(fact(ii - 1)*2**(ii - 1))))* &
                        &((m*omega)/(pi*hbar))* &
                        &exp(-(m*omega*grid_points**2)/(2*hbar))* &
                        &evalHermitePoly(sqrt((m*omega)/hbar)*grid_points, ii - 1)
        eigenf(:, ii) = eigenf(:, ii) / norm2(real(eigenf(:, ii)))
    end do

    allocate (prob(N, k))

    call ComputeProb(eigenf, prob)

    eigenf_filename = "ef_t_"//trim(str_i(k))//"_"//trim(str_i(L))//"_"//trim(str_r_e(dx))//"_"//trim(str_r_e(omega))//"_"&
                      &//trim(str_r_d(m))//"_"//trim(str_r_d(hbar))//".txt"
    prob_filename = "pr_t_"//trim(str_i(k))//"_"//trim(str_i(L))//"_"//trim(str_r_e(dx))//"_"//trim(str_r_e(omega))//"_"&
                    &//trim(str_r_d(m))//"_"//trim(str_r_d(hbar))//".txt"

    call WriteEigenvectors(real(eigenf), grid_points, eigenf_filename)
    call WriteEigenvectors(prob, grid_points, prob_filename)

end program hermite_polynomals
