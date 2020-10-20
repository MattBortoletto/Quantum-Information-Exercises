module matrices

    implicit none

    type dmatrix
        integer, dimension(2) :: dim 
        double complex, dimension(:, :), allocatable :: elem
        double complex :: tr, det 
    end type

    interface operator (.init.)
        module procedure MatInit
    end interface operator (.init.)

    interface operator (.trace.)
        module procedure MatTrace
    end interface operator (.trace.)

    interface operator (.adj.)
        module procedure MatAdj
    end interface operator (.adj.) 

    contains

    ! function to initialize a matrix 
    function MatInit(dimensions) result(M)

        integer, dimension(2), intent(in) :: dimensions 
        type(dmatrix) :: M 
        if ((dimensions(1) .le. 1) .and. (dimensions(2) .le. 1)) then 
            print *, "[MatInit warning] Non valid matrix dimension."
        else
            allocate(M%elem(dimensions(1), dimensions(2)))
            M%dim = dimensions
            M%elem = 0.d0
            M%tr = 0.d0
            M%det = 0.d0
        end if  

        return 
        
    end function MatInit

    ! function that computes the trace of a matrix
    function MatTrace(M) result(trace) 

        type(dmatrix), intent(in) :: M 
        double complex :: trace 
        integer :: ii 

        trace = 0.d0

        if (M%dim(1) .ne. M%dim(2)) then    ! check if the matrix M is a square matrix
            print *, "[MatTrace warning] The matrix is not square. The trace will be set to zero." 
            trace = 0.d0
        else 
            do ii = 1, M%dim(1) 
                trace = trace + M%elem(ii, ii) 
            end do 
        end if 

        return 

    end function MatTrace

    ! function which computes the adjoint matrix and its trace, determinant and dimension
    function MatAdj(M) result(Madj)

        type(dmatrix), intent(in) :: M 
        type(dmatrix) :: Madj 

        Madj%tr = conjg(M%tr) 
        Madj%det = conjg(M%det) 
        Madj%dim = M%dim(2:1:-1)
        Madj%elem = transpose(conjg(M%elem))

        return 

    end function MatAdj

    ! subroutine that writes the results in a text file
    subroutine WriteResults(M, filename) 

        type(dmatrix) :: M
        character*8 :: filename 
        integer :: ll

        open(unit = 2, file = filename, action = "write", status = "replace") 
        
        write(2, *) "MATRIX:"
        do ll = 1, M%dim(1)
            write(2, *) M%elem(ll, :) 
        end do
        write(2, *) " "
        write(2, *) "DIMENSION: rows =",  M%dim(1), ",  columns = ", M%dim(2)
        write(2, *) " "
        write(2, *) "TRACE:", M%tr
        write(2, *) " "
        write(2, *) "DETERMINANT:", M%det

        close(2)

        return 

    end subroutine WriteResults

end module matrices


program main

    use matrices 

    implicit none

    ! A: original matrix
    ! Aadj: adjoint matrix
    type(dmatrix) :: A, Aadj
    integer :: nrows, ncols
    integer, dimension(2) :: d 
    real, dimension(:, :), allocatable :: RealPart, ImPart
    character*8 :: file1, file2
    character*3 :: answ

    nrows = 0

    ! enter the dimension of the matrix
    print *, "Please enter the dimension of the matrix [nrows, ncols]:"
    read (*, *) nrows, ncols
    
    ! ask to enter the dimension until nrows and ncols are greater than 1
    do while ((nrows .le. 1) .or. (ncols .le. 1))
        print *, "Non valid dimension! The number of rows and columns must be greater than 1."
        print *, "Please enter the dimension of the matrix [nrows, ncols]:"
        read (*, *) nrows, ncols
        ! in case of non-square matrix ask if it is ok 
        if (nrows .ne. ncols) then 
            print *, "Warining! The matrix is not square, so the trace is not defined and will be set to zero."
            print *, "Do you want to continue? [y/n]" 
            read *, answ 
            ! if the user answers 'y' then in the MatTrace function the trace will be set to zero
            if (answ == "y") then 
                continue
            ! if the user answers 'n' then he will be asked to enter the dimension of a square matrix
            else if (answ == "n") then 
                do while (nrows .ne. ncols) 
                    print *, "Please enter the dimension of a square matrix [nrows, ncols]:"
                    read (*, *) nrows, ncols
                    if (nrows .ne. ncols) then 
                        print *, "The matrix is not square."
                    end if 
                end do
            else
                print *, "Non valid answer. The program will stop."
                stop 
            end if 
        end if 
    end do

    allocate(RealPart(nrows, ncols))
    allocate(ImPart(nrows, ncols))

    ! fill the matrices with random numbers between 0 and 1
    call random_number(RealPart)
    call random_number(ImPart)

    d(1) = nrows
    d(2) = ncols

    A = .init.(d)
    A%elem = cmplx(RealPart, ImPart)
    A%tr = .trace.(A)
    Aadj = .adj.(A)

    ! save the results
    file1 = "M.txt"
    file2 = "Madj.txt"
    call WriteResults(A, file1)
    call WriteResults(Aadj, file2)

    deallocate(A%elem, Aadj%elem)

    stop 
    
end program main