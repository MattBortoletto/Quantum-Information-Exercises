module variables
    integer*2 :: n1, n2
    integer*4 :: n3, n4
    real*4 :: n5, n6, pi4
    real*8 :: n7, n8, pi8
end module variables 


program MySum

    use variables
    implicit none

    !n1 = 2000000
    !n2 = 1
    print *, "With integer*2 2.000.000 + 1 leads to an overflow error"

    n3 = 2000000
    n4 = 1
    print *, "With integer*4, 2.000.000 + 1 = ", n3 + n4

    pi4 = acos(-1.)
    n5 = pi4 * 10e32
    n6 = sqrt(2.) * 10e21
    print *, "In single precision pi*10e32 + sqrt(2)*10e21 =", n5 + n6

    pi8 = 4.D0 * datan(1.D0)
    n7 = pi8 * 10e32
    n8 = sqrt(2.) * 10e21
    print *, "In single precision pi*10e32 + sqrt(2)*10e21 =", n7 + n8

    stop

end program MySum