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
        !real*4, dimension(2) :: GaussReIm
        ! flag to choose the type of matrix to initialize
        character(1), intent(in) :: which_matrix

        if (which_matrix == "h") then
            do ii = 1, d
                ! first, fill the diagonal
                ! in order to be hermitian, the matrix must have its
                ! diagonal elements with no imaginary part
                call random_number(RealPart)
                !call random_number(ImPart)
                !GaussReIm = BoxMuller(RealPart, ImPart)
                !matr(ii, ii) = complex(GaussReIm(1), 0)
                RealPart = 2.0*RealPart - 1.0
                matr(ii, ii) = complex(RealPart, 0)
                ! then fill the rest
                ! in order to be hermitian, the matrix mush be such
                ! that matr(ii, jj) = conj(matr(ii, jj)) 
                do jj = ii + 1, d 
                    call random_number(RealPart)
                    call random_number(ImPart)
                    !GaussReIm = BoxMuller(RealPart, ImPart)
                    !matr(ii, jj) = complex(GaussReIm(1), GaussReIm(2))
                    RealPart = 2.0*RealPart - 1.0
                    ImPart = 2.0*ImPart - 1.0
                    matr(ii, jj) = complex(RealPart, ImPart)
                    matr(jj, ii) = conjg(matr(ii, jj))
                end do
            end do 
        ! if we choose "d", the matrix will be real and diagonal 
        else if (which_matrix == "d") then 
            matr = complex(0.d0, 0.d0)
            do ii = 1, d
                call random_number(RealPart)
                !call random_number(ImPart)
                !GaussReIm = BoxMuller(RealPart, ImPart)
                !matr(ii, ii) = GaussReIm(1)
                RealPart = 2.0*RealPart - 1.0
                matr(ii, ii) = RealPart
            end do
        else 
            print *, "Non-valid input!" 
            stop 
        end if 

        return 

    end function MatrixInit


    subroutine ComputeEigenvalues(matr, eig, info)

        ! this subroutine computes the eigenvalues of a complex hermitian matrix
        ! using the LAPACK subroutine 'cheev', which takes this arguments:
        ! - jobz = "N": compute eigenvalues only
        ! - uplo = "U": upper triangle of A is stored
        ! - N = size(matr, 1): order of the matrix
        ! - a = matr: input matrix
        ! - lda = size(matr, 1): leading dimension of the matrix
        ! - w: if info = 0, the eigenvalues in ascending order
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


    function ComputeSpacings(eig) result(spacings)

        ! vector to store the eigenvalues
        real*4, dimension(:), intent(in) :: eig 
        ! vectors to store the spacings and the normalized spacings
        real*4, dimension(size(eig, 1)-1) :: spacings
        ! variable to loop
        integer :: ii 

        ! compute the spacings
        do ii = 1, size(eig, 1) - 1
            spacings(ii) = eig(ii+1) - eig(ii)
        end do

        return 

    end function ComputeSpacings


    function ComputeNormSpacings(eig) result(norm_spacings)

        ! this subroutine computes the normalized spacings between eigenvalues.
        ! in order to do that we need to compute the difference between adjacent
        ! eigenvalues and the average of these spacings

        ! vector to store the eigenvalues
        real*4, dimension(:), intent(in) :: eig 
        ! vectors to store the spacings and the normalized spacings
        real*4, dimension(size(eig, 1)-1) :: spacings, norm_spacings
        ! mean of the spicings
        real*4 :: mean 
        
        ! compute the spacings
        spacings = ComputeSpacings(eig)
        
        ! compute the mean of the spacings 
        mean = sum(spacings) / size(spacings, 1)

        ! normalize the spacings
        norm_spacings = spacings / mean 

        return 

    end function ComputeNormSpacings


    function ComputeNormSpacingsLocal(eig, range) result(spacings_local)

        ! this subroutine computes the normalized spacings between eigenvalues.
        ! in order to do that we need to compute the difference between adjacent
        ! eigenvalues and the average of these spacings

        ! vector which in our case stores the normalized spacings 
        real*4, dimension(:) :: eig
        real*4, dimension(size(eig)-1) :: spacings, local_avg, spacings_local 
        ! variable to loop
        integer :: jj
        ! window size 
        integer, intent(in) :: range 

        if (range .gt. size(eig)) then 
            print *, "The window size is too big!"
            stop 
        end if 

        ! compute the spacings 
        spacings = ComputeSpacings(eig)
        
        ! compute the local averages
        do jj = 1, size(spacings, 1) 
            local_avg(jj) = sum(spacings( max(1, jj-range):min(jj+range, size(spacings)) ))
            local_avg(jj) = local_avg(jj) / (min(jj+range, size(spacings)) - max(1, jj-range))
        end do 

        ! normalize 
        spacings_local = spacings / local_avg

        ! print *, spacings_local
        ! print *, " "

        return  

    end function ComputeNormSpacingsLocal


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
        real*4, dimension(:), allocatable :: dist 

        allocate(right_edge(nbins))
        allocate(counts(nbins))

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
        dist = real(counts) / (sum(counts)*bin_size)

        deallocate(right_edge)
        deallocate(counts)

    end subroutine ComputePDF


    function ComputeR(spacings) result(r_mean) 

        real*4, dimension(:) :: spacings
        real*4, dimension(size(spacings-1)) :: r 
        real*4 :: r_mean 
        integer :: ii 

        do ii = 1, size(r)
            r(ii) = min(spacings(ii), spacings(ii+1)) / max(spacings(ii), spacings(ii+1))
        end do 

        r_mean = sum(r) / size(r) 

        return

    end function ComputeR


    function str(k) result(k_str)

        ! convert an integer to string
        integer, intent(in) :: k
        character(30) :: k_str

        write (k_str, *) k
        k_str = adjustl(k_str)

    end function str

end module rand_matrix



program RandomMatrix 
    
    use rand_matrix

    implicit none

    ! N: dimension of the matrix
    ! ii, jj, trial: variables to cycle
    ! ntrials: number of matrices generated to compute the distribution 
    !          distribution of the spacings 
    ! n_bins: number of bins for the histogram
    ! div: number of local spacings groups -------------------------------------- TOGLIERE
    integer :: N, ii, jj, trial, ntrials, info, n_bins !, div
    ! matrix
    complex, dimension(:,:), allocatable :: M
    ! eig: array to store the eigenvalues
    ! norm_spacings: array to store the normalized spacings of a 
    !                single matrix 
    ! all_s: vector to store all the normalized spacings 
    real*4, dimension(:), allocatable :: eig, norm_spacings, all_s
    ! which_matrix: flag to choose the type of matrix (hermitian or real diagonal) 
    ! which_spacing: flag to choose the type of spacings (normalized
    !                with respect to the global or local mean)
    character(1) :: which_matrix, which_spacing
    ! hist_bins: array to store the bin centers of the distribution
    ! distribution: array to store the distribution points 
    real*4, dimension(:), allocatable :: hist_bins, distribution 
    ! filename: name for the text file in which the results are stored 
    character(30) :: filename
    ! variable to store <r>
    real*4 :: r_mean 
    ! div
    integer, dimension(5) :: div
    ! local_norm_spacings: array to store the locally normalized spacings
    !                      of a single matrix for different locality levels
    real*4, dimension(:, :), allocatable :: local_norm_spacings
    real*4, dimension(:, :), allocatable :: all_local_s

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

    ! ask the user to choose between global or local mean  
    print *, "Do you want to consider the global average or local averages? [g/l]"
    read *, which_spacing
    if ((which_spacing .ne. "g") .and. (which_spacing .ne. "l")) then 
        print *, "Invalid input."
        stop
    end if 
    ! if (which_spacing == "l") then
    !     print *, "Please enter the number of divisions D to form the local spacings (N/D):"
    !     read *, div 
    ! end if 

    div = (/100, 50, 10, 5, 1/)


    ! set the number random matrices which are used to
    ! compute the spacings distribution 
    print *, "How many matrices do you want to generate?"
    read *, ntrials 

    ! allocate the memory
    allocate(M(N, N))
    allocate(eig(N))
    allocate(norm_spacings(N-1)) ! HO MESSO -1
    allocate(all_s((N-2)*ntrials))
    allocate(local_norm_spacings(N-1, size(div))) ! HO MESSO -1
    allocate(all_local_s((N-2)*ntrials, size(div)))

    do trial = 1, ntrials

        print *, "Computing the spacings for matrix number", trial

        ! use the info flag of the subroutine 'cheev' to check if the
        ! diagonalization has been successful
        info = 1
        do while (info .ne. 0)
            ! initialize the matrix 
            M = MatrixInit(N, which_matrix)
            ! call the subroutine to compute the eigenvalues 
            call ComputeEigenvalues(M, eig, info)
        end do 

        ! call the function to compute the spacings
        if (which_spacing == "g") then  
            norm_spacings = ComputeNormSpacings(eig)
            ! add the spacings we just computed to the vector of all 
            ! the spacings 
            all_s(1+(trial-1)*(N-2):trial*(N-2)) = norm_spacings
        else if (which_spacing == "l") then 
            do ii = 1, size(div)
                local_norm_spacings(:, ii) = ComputeNormSpacingsLocal(eig, N/div(ii))
                all_local_s(1+(trial-1)*(N-2):trial*(N-2), ii) = local_norm_spacings(:, ii)
            end do 
        end if 

        ! add the spacings we just computed to the vector of all 
        ! the spacings 
        ! all_s(1+(trial-1)*(N-2):trial*(N-2)) = norm_spacings

    end do 

    ! compute the <r>
    if (which_spacing == "g") then 
        r_mean = ComputeR(all_s)
        print *, "<r> =", r_mean
    else if (which_spacing == "l") then 
        do ii = 1, size(div)
            print *, "Locality level = N /", div(ii)
            r_mean = ComputeR(local_norm_spacings(:, ii))
            print *, "<r> =", r_mean
        end do 
    end if   

    ! use Rice Rule to compute the optimal number of bins
    ! nbins = 2*N^{1/3}
    n_bins = int(2.0*(ntrials*(N-2))**(1.0/3.0))

    allocate(hist_bins(n_bins))
    allocate(distribution(n_bins))

    ! call the subroutine which computes the pdf
    ! call ComputePDF(all_s, n_bins, distribution, hist_bins)

    ! save the results in a text file according to the chosen options
    ! if ((which_matrix == "h") .and. (which_spacing == "g")) then 
    !     filename = "herm_spacings.txt"
    !     open(10, file=filename, status='replace')
    !     do ii = 1, size(hist_bins)
    !         write(10, *) hist_bins(ii), distribution(ii)
    !     end do
    !     write(10, *)
    !     close(10) 
    ! else if ((which_matrix == "d") .and. (which_spacing == "g")) then
    !     filename = "diag_spacings.txt"
    !     open(10, file=filename, status='replace')
    !     do ii = 1, size(hist_bins)
    !         write(10, *) hist_bins(ii), distribution(ii)
    !     end do
    !     write(10, *)
    !     close(10)
    ! end if 

    ! save the results in a text file 
    ! choose the name of the file according to the options the user chose
    if ((which_matrix == "h") .and. (which_spacing == "g")) then 
        call ComputePDF(all_s, n_bins, distribution, hist_bins)
        filename = "herm_spacings.txt"
        open(10, file=filename, status='replace')
        do ii = 1, size(hist_bins)
            write(10, *) hist_bins(ii), distribution(ii)
        end do
        write(10, *)
        close(10) 
    else if ((which_matrix == "h") .and. (which_spacing == "l")) then 
        do jj = 1, size(div)
            call ComputePDF(all_local_s(:, jj), n_bins, distribution, hist_bins)
            filename = "herm_loc_aver_"//trim(str(div(jj)))//".txt"
            open(10, file=filename, status='replace')
            do ii = 1, size(hist_bins)
                write(10, *) hist_bins(ii), distribution(ii)
            end do
            write(10, *)
            close(10)
        end do 
    else if ((which_matrix == "d") .and. (which_spacing == "g")) then
        call ComputePDF(all_s, n_bins, distribution, hist_bins)
        filename = "diag_spacings.txt"
        open(10, file=filename, status='replace')
        do ii = 1, size(hist_bins)
            write(10, *) hist_bins(ii), distribution(ii)
        end do
        write(10, *)
        close(10)
    else if ((which_matrix == "d") .and. (which_spacing == "l")) then
        do jj = 1, size(div)
            call ComputePDF(all_local_s(:, jj), n_bins, distribution, hist_bins)
            filename = "diag_loc_aver_"//trim(str(div(jj)))//".txt"
            open(10, file=filename, status='replace')
            do ii = 1, size(hist_bins)
                write(10, *) hist_bins(ii), distribution(ii)
            end do
            write(10, *)
            close(10)
        end do 
    end if 

    ! open(10, file=filename, status='replace')
    ! do ii = 1, size(hist_bins)
    !     write(10, *) hist_bins(ii), distribution(ii)
    ! end do
    ! write(10, *)
    ! close(10) 
    
    deallocate(m)
    deallocate(eig)
    deallocate(distribution)
    deallocate(hist_bins)
    
end program RandomMatrix
