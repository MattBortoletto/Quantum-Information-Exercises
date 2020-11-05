module rand_matrix

    implicit none 

contains

    function BoxMuller(a, b) result(x) 

        ! this function takes in input two values 'a' and 'b'
        ! and computes a pair of gaussian distributed variables
        ! with zero mean and unitary variance using the Box-Muller
        ! algorithm  
        
        real*4, intent(in) :: a, b 
        real*4, dimension(2) :: x 
        real*4 :: TwoPi 

        TwoPi = 8.0d0*atan(1.0)

        x(1) = sqrt(-2.0*log(a))*cos(TwoPi*b)
        x(2) = sqrt(-2.0*log(a))*cos(TwoPi*b)
        
        return 

    end function BoxMuller


    function MatrixInit(d, which_matrix) result(matr) 

        ! dimension of the matrix
        integer :: d  
        ! matrix elements
        complex, dimension(d, d) :: matr
        ! variables to loop
        integer :: ii, jj 
        ! real part and imaginary part
        real*4 :: RealPart, ImPart
        ! gaussian distributed real and imaginary parts
        real*4, dimension(2) :: GaussReIm
        ! flag to choose the type of matrix to initialize
        character(1), intent(in) :: which_matrix

        if (which_matrix == "h") then
            do ii = 1, d
                ! first, fill the diagonal
                ! in order to be hermitian, the matrix must have its
                ! diagonal elements with no imaginary part
                call random_number(RealPart)
                call random_number(ImPart)
                GaussReIm = BoxMuller(RealPart, ImPart)
                matr(ii, ii) = complex(GaussReIm(1), 0)
                ! then fill the rest
                ! in order to be hermitian, the matrix mush be such
                ! that matr(ii, jj) = conj(matr(ii, jj)) 
                do jj = ii + 1, d 
                    call random_number(RealPart)
                    call random_number(ImPart)
                    GaussReIm = BoxMuller(RealPart, ImPart)
                    matr(ii, jj) = complex(GaussReIm(1), GaussReIm(2))
                    matr(jj, ii) = conjg(matr(ii, jj))
                end do
            end do 
        ! if we choose "d", the matrix will be real and diagonal 
        else if (which_matrix == "d") then 
            matr = complex(0.d0, 0.d0)
            do ii = 1, d
                call random_number(RealPart)
                call random_number(ImPart)
                GaussReIm = BoxMuller(RealPart, ImPart)
                matr(ii, ii) = GaussReIm(1)
            end do
        else 
            print *, "Non-valid input!" 
            stop 
        end if 

        return 

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
        real*4, dimension(size(matr, 1)) :: eig 
        real*4, dimension(:), allocatable :: rwork
        complex, dimension(:), allocatable :: work
        character(1) :: jobz, uplo
        integer :: n, lda, lwork, info
        
        n = size(matr, 1)
        lda = size(matr, 1)
        lwork = 2*size(matr, 1) - 1
        jobz = "N"
        uplo = "U"
        allocate(rwork(3*size(matr, 1)-2))
        allocate(work(2*size(matr, 1)-1))
        
        call cheev(jobz, uplo, n, matr, lda, eig, work, lwork, rwork, info)

        deallocate(rwork, work)

        return 

    end subroutine ComputeEigenvalues


    function ComputeSpacings(eig) result(norm_spacings)

        ! this subroutine computes the normalized spacings between eigenvalues.
        ! in order to do that we need to compute the difference between adjacent
        ! eigenvalues and the average of these spacings

        ! vector to store the eigenvalues
        real*4, dimension(:), intent(in) :: eig 
        ! vectors to store the spacings and the normalized spacings
        real*4, dimension(size(eig, 1)-1) :: spacings, norm_spacings
        ! mean of the spicings
        real*4 :: mean 
        ! variable to loop
        integer :: ii
        
        do ii = 1, size(eig, 1) - 1
            spacings(ii) = eig(ii+1) - eig(ii)
        end do
        
        mean = sum(spacings) / size(spacings, 1)

        norm_spacings = spacings / mean 

        return 

    end function ComputeSpacings


    function ComputeSpacingsLocal(eig, range) result(norm_spacings_local)

        ! this subroutine computes the normalized spacings between eigenvalues.
        ! in order to do that we need to compute the difference between adjacent
        ! eigenvalues and the average of these spacings

        ! vector to store the eigenvalues
        real*4, dimension(:), intent(in) :: eig 
        ! vectors to store the local spacings
        real*4, dimension(size(eig, 1)-1) :: spacings
        ! vectors to store the local normalized spacings
        real*4, dimension(:), allocatable :: norm_spacings_local
        ! variables to loop
        integer :: ii, jj
        ! window size 
        integer, intent(in) :: range 
        
        do ii = 1, size(eig, 1) - 1
            spacings(ii) = eig(ii+1) - eig(ii)
        end do 

        allocate(norm_spacings_local(size(spacings) - range + 1))

        do jj = 1, size(spacings) - range + 1
            norm_spacings_local(jj) = sum(spacings(jj:jj+range-1)) 
            norm_spacings_local(jj) = norm_spacings_local(jj) / range 
        end do 

        return  

    end function ComputeSpacingsLocal


    subroutine ComputePDF(x, nbins, dist, bin_centers)

        ! input vector 
        real*4, dimension(:), intent(in) :: x 
        ! number of bins 
        integer, intent(in) :: nbins
        ! bin width and temporary variable which is used to compute  
        ! the right edges of the bins
        real*4 :: bin_size, bin_increment   
        ! vectors to store the right edges of the bins and the bin centers
        real*4, dimension(:), allocatable :: right_edge, bin_centers
        ! variables to loop 
        integer :: ii, jj  
        ! vector to store the counts for the histogram
        integer, dimension(:), allocatable :: counts
        ! vectors which store the normalized counts and the pdf 
        real*4, dimension(:), allocatable :: norm_counts, dist 

        allocate(right_edge(nbins))
        allocate(counts(nbins))
        allocate(norm_counts(nbins))

        ! compute the size of the bins 
        bin_size = (maxval(x) - minval(x)) / nbins
        
        ! compute the right edges of the bins
        bin_increment = minval(x)
        do ii = 1, nbins
            bin_increment = bin_increment + bin_size
            right_edge(ii) = bin_increment
        end do

        ! compute the bin centers
        bin_centers = right_edge - (bin_size/2)

        ! fill the histogram
        counts = 0
        do ii = 1, size(x, 1)                       
            do jj = 1, nbins - 1                
                if (x(ii) .le. right_edge(jj)) then
                    counts(jj) = counts(jj) + 1
                    exit
                end if           
            end do                         
            if (x(ii) .ge. right_edge(nbins-1)) counts(nbins) = counts(nbins) + 1
        end do

        ! normalization 
        norm_counts = real(counts) / sum(counts)

        ! compute the pdf 
        do ii=1, nbins 
            dist(ii) = norm_counts(ii) / (bin_size)
        end do

        deallocate(right_edge)
        deallocate(counts)
        deallocate(norm_counts) 

    end subroutine ComputePDF

end module rand_matrix



program RandomMatrix 
    
    use rand_matrix

    implicit none

    ! dimension of the matrix (N) and variables to cycle
    integer :: N, ii, trial, ntrials
    ! matrix
    complex, dimension(:,:), allocatable :: M
    ! array to store the eigenvalues, the normalized spacings of a 
    ! single matrix and all the normalized spacings 
    real*4, dimension(:), allocatable :: eig, norm_spacings, all_s
    ! flag to choose the type of matrix (hermitian or real diagonal) 
    character(1) :: which_matrix, which_spacing
    ! histogram bin centers and counts 
    real*4, dimension(:), allocatable :: hist_bins, hist_counts
    ! number of bins for the histogram
    integer :: n_bins
    ! name for the text file in which the results are stored 
    character(21) :: filename 

    ! ask the user to enter the dimension of the matrix
    print *, "Please enter the dimension of the matrix: "
    read *, N 

    ! ask the user to choose the type of matrix
    print *, "Do you want a hermitian or a diagonal matrix? [h/d]"
    read *, which_matrix
    if ((which_matrix .ne. "h") .and. (which_matrix .ne. "d")) then 
        print *, "Invalid input."
        stop
    end if 

    ! ask the user to choose between global or local spacings 
    print *, "Do you want to use global or local spacings? [g/l]"
    read *, which_spacing
    if ((which_spacing .ne. "g") .and. (which_spacing .ne. "l")) then 
        print *, "Invalid input."
        stop
    end if 

    ! allocate the memory
    allocate(M(N, N))
    allocate(eig(N))
    allocate(norm_spacings(N))
    allocate(all_s((N-2)*100))

    ! set the number random matrices which are used to
    ! compute the spacings distribution 
    ntrials = 100

    do trial = 1, ntrials

        ! initialize the matrix 
        M = MatrixInit(N, which_matrix)

        ! call the subroutine to compute the eigenvalues 
        call ComputeEigenvalues(M, eig)

        ! call the function to compute the spacings 
        if (which_spacing == "g") then 
            norm_spacings = ComputeSpacings(eig)
        else if (which_spacing == "l") then 
            norm_spacings = ComputeSpacingsLocal(eig, 5)
        end if 

        ! add the spacings we just computed to the vector of all 
        ! the spacings 
        all_s(1+(trial-1)*(N-2):trial*(N-2)) = norm_spacings

        print *, "Computing the spacings for matrix number", trial

    end do 

    ! use Rice Rule to compute the optimal number of bins
    ! nbins = 2*N^{1/3}
    n_bins = int(2.0*(ntrials*(N-2))**(1.0/3.0))
    print *, n_bins

    allocate(hist_bins(n_bins))
    allocate(hist_counts(n_bins))

    ! call the subroutine which computes the pdf
    call ComputePDF(all_S, n_bins, hist_counts, hist_bins)

    ! save the results in a text file 
    if ((which_matrix == "h") .and. (which_spacing == "g")) then 
        filename = "herm_glo_spacings.txt"
    else if ((which_matrix == "h") .and. (which_spacing == "l")) then 
        filename = "herm_loc_spacings.txt"
    else if ((which_matrix == "d") .and. (which_spacing == "g")) then
        filename = "diag_glo_spacings.txt"
    else if ((which_matrix == "d") .and. (which_spacing == "l")) then
        filename = "diag_loc_spacings.txt"
    end if 
    open(10, file=filename, status='replace')
    do ii = 1, size(hist_bins)
        write(10, *) hist_bins(ii), hist_counts(ii)
    end do
    write(10, *)
    
    close(10) 
    
    deallocate(m)
    deallocate(eig)
    deallocate(hist_counts)
    deallocate(hist_bins)
    
end program RandomMatrix
