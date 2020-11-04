module herm_rand_matrix

    implicit none 

contains

    function MatrixInit(d, which_matrix) result(matr) 

        ! dimension of the matrix
        integer :: d  
        ! matrix elements
        complex, dimension(d, d) :: matr
        ! variables to loop
        integer :: ii, jj 
        ! real part and imaginary part
        real :: RealPart, ImPart, tmp
        ! flag
        character(1), intent(in) :: which_matrix

        if (which_matrix == "h") then
            do ii = 1, d
                ! first, fill the diagonal
                ! in order to be hermitian, the matrix must have its
                ! diagonal elements with no imaginary part
                call random_number(RealPart)
                matr(ii, ii) = complex(RealPart, 0)
                ! then fill the rest
                ! in order to be hermitian, the matrix mush be such
                ! that matr(ii, jj) = conj(matr(ii, jj)) 
                do jj = ii + 1, d 
                    call random_number(RealPart)
                    call random_number(ImPart)
                    matr(ii, jj) = complex(RealPart, ImPart)
                    matr(jj, ii) = conjg(matr(ii, jj))
                end do
            end do 
        else if (which_matrix == "d") then 
            matr = complex(0.d0, 0.d0)
            do ii = 1, d
                call random_number(tmp)
                matr(ii, ii) = tmp
            end do
        else 
            print *, "Non-valid input!" 
            stop 
        end if 

    end function MatrixInit


    subroutine ComputeEigenvalues(matr, eig)

        ! this subroutine computes the eigenvalues of a complex hermitian matrix
        ! using the LAPACK subroutine 'cheev', which takes this arguments:
        ! - jobz = "N": compute eigenvalues only
        ! - uplo = "U": upper triangle of A is stored
        ! - N = size(matr, 1): order of the matrix
        ! - a = matr: input matrix
        ! - lda = size(matr, 1): leading dimension of the matrix
        ! - w = : if info = 0, the eigenvalues in ascending order
        ! - work: complex array of dimension (max(1,lwork))
        !         on exit, if info = 0, work(1) returns the optimal lwork
        ! - lwork: length of the array work
        ! - rwork: workspace
        !          real array, dimension (max(1, 3*n-2))
        ! - info: output, if info = 0 the exit is successful

        complex, dimension(:, :), intent(in) :: matr
        real, dimension(size(matr, 1)) :: eig 
        real, dimension(:), allocatable :: rwork
        complex, dimension(:), allocatable :: work
        character(1) :: jobz, uplo
        integer :: n, lda, lwork, info
        
        n = size(matr, 1)
        lda = size(matr, 1)
        lwork = 2*size(matr, 1)-1
        jobz = "N"
        uplo = "U"
        allocate(rwork(3*size(matr, 1)-2))
        allocate(work(2*size(matr, 1) - 1))
        
        call cheev(jobz,uplo,n,matr,lda,eig,work,lwork,rwork,info)

        deallocate(rwork, work)

    end subroutine ComputeEigenvalues


    function ComputeSpacings(eig) result(norm_spacings)

        ! this subroutine computes the normalized spacings between eigenvalues.
        ! in order to do that we need to compute the difference between adjacent
        ! eigenvalues and the average of these spacings
        real, dimension(:), intent(in) :: eig 
        real, dimension(size(eig, 1)-1) :: spacings, norm_spacings
        real :: mean 
        integer :: ii
        
        ! print *, "eig", eig
        
        do ii = 1, size(eig, 1) - 1
            spacings(ii) = eig(ii+1) - eig(ii)
        end do

        ! print *, "spacings", spacings
        
        mean = sum(spacings) / size(spacings, 1)

        ! print *, "mean", mean 

        norm_spacings = spacings / mean 

        ! print *, "norm spacings", norm_spacings

    end function ComputeSpacings


    function ComputeSpacingsLocal(eig, range) result(norm_spacings_local)

        ! this subroutine computes the normalized spacings between eigenvalues.
        ! in order to do that we need to compute the difference between adjacent
        ! eigenvalues and the average of these spacings
        real, dimension(:), intent(in) :: eig 
        real, dimension(size(eig, 1)-1) :: spacings
        real, dimension(:), allocatable :: norm_spacings_local
        integer :: ii, jj
        integer, intent(in) :: range 
        
        do ii = 1, size(eig, 1) - 1
            spacings(ii) = eig(ii+1) - eig(ii)
        end do 

        allocate(norm_spacings_local(size(spacings) - range + 1))

        do jj = 1, size(spacings) - range + 1
            norm_spacings_local(jj) = sum(spacings(jj: jj+range-1)) 
            norm_spacings_local(jj) = norm_spacings_local(jj) / range 
        end do 

        !print *, "norm spacings local after normalization", norm_spacings_local

    end function ComputeSpacingsLocal


    subroutine hist(x, nbins, counts, bin_centers)

        ! input vector 
        real, dimension(:), intent(in) :: x 
        ! number of bins 
        integer, intent(in) :: nbins
        real :: bin_size, bin_increment    
        real, dimension(:), allocatable :: right_edge
        real, dimension(:), allocatable :: bin_centers      
        integer :: ii, jj  
        integer, dimension(:), allocatable :: counts

        allocate(right_edge(nbins))
        allocate(bin_centers(nbins))
        allocate(counts(nbins))

        bin_size = (maxval(x) - minval(x)) / nbins
        ! print *, "bin size", bin_size
        
        bin_increment = minval(x)
        do ii = 1, nbins
            bin_increment = bin_increment + bin_size
            right_edge(ii) = bin_increment
        end do

        bin_centers = right_edge - (bin_size/2)

        ! print *, "right edges", right_edge 

        counts = 0
        do ii = 1, size(x,1)                       
            do jj = 1, nbins-1
                ! print *, "x_i", x(ii)                  
                if (x(ii) .le. right_edge(jj)) then
                counts(jj) = counts(jj) + 1
                exit
                end if           
            end do                         
            if (x(ii) >= right_edge(nbins-1)) counts(nbins) = counts(nbins) + 1
        end do

        ! print *, "from the inside ... ", counts

        return 

    end subroutine hist

end module herm_rand_matrix



program eigenproblem 
    
    use herm_rand_matrix

    implicit none

    ! dimension of the matrix
    integer :: N
    ! matrix
    complex, dimension(:,:), allocatable :: M
    ! array to store the eigenvalues 
    real, dimension(:), allocatable :: eig, norm_spacings, norm_spacings_local
    ! flag 
    character(1) :: which_matrix

    real, dimension(:), allocatable :: hist_bins 
    integer, dimension(:), allocatable :: hist_counts

    integer :: n_bins

    ! ask the user to enter the dimension of the matrix
    print *, "Please enter the dimension of the matrix: "
    read *, N 

    print *, "Do you want a hermitian or a diagonal matrix?"
    read *, which_matrix

    ! allocate the memory
    allocate(M(N, N))
    allocate(eig(N))
    allocate(norm_spacings(N))
    ! allocate(hist_bins(n_bins))
    ! allocate(hist_counts(n_bins))

    ! initialize the matrix 
    M = MatrixInit(N, which_matrix)

    ! print *, M

    ! call the subroutine to compute the eigenvalues 
    call ComputeEigenvalues(M, eig)

    ! call the function to compute the spacings 
    norm_spacings = ComputeSpacings(eig)

    print *, "norm spacings", norm_spacings

    norm_spacings = norm_spacings(1:size(norm_spacings, 1) - 1)

    print *, "norm spacings", norm_spacings

    ! call the function to compute the local spacings 
    norm_spacings_local = ComputeSpacingsLocal(eig, 5)

    ! save the results in a text file 
    ! open(10, file='results.txt', status='replace')
    ! do ii = 1, n
    !     write(10, *) ii, ????????????????????(ii)
    !     ! 10 format(I3, '   ', f14.8)
    ! end do
    ! write(10, *)
    
    ! close(10) 

    n_bins = 32

    call hist(norm_spacings, n_bins, hist_counts, hist_bins)

    print *, "counts", hist_counts
    print *, "bins", hist_bins 
    
    deallocate(m)
    deallocate(eig)
    ! deallocate(hist_counts)
    ! deallocate(hist_bins)
    
end program eigenproblem
