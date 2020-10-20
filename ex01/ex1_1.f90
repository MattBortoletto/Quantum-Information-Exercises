module variables1

    integer*2 :: a, b 

end module variables1


program MyFirstProgram

    use variables1

    implicit none

    a = 15
    b = 18

    print *, "Hello World! 15 + 18 =", a + b 

    stop
    
end program MyFirstProgram