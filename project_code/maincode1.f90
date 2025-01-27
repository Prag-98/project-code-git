program main1 !CODE USING MODULES - GENERATING RANDOM MATRIX & DIAGONALIZING IT
    use, intrinsic :: iso_fortran_env, only: a=>real64
    use mymodule1
    implicit none

    integer :: i, io1, io2
    integer :: lwork, liwork, lrwork, info
    integer, parameter :: n=12000, nmax=(1+5*n+2*n**2)
    real(a), parameter :: error=1e-8

    complex(a), allocatable :: mat(:,:), hermmat(:,:), egnvecs(:,:)
    integer, allocatable :: iwork(:)
    real(a), allocatable :: egnvls(:), rwork(:)
    complex(a), allocatable :: work(:)
    complex(a), allocatable :: mat1(:,:), mat2(:,:), matcheck(:,:), hermmat1(:,:), Q(:,:), diff(:,:), egndiag(:,:)

    call random_mtx(n,mat)
    mat = mat*100

    hermmat = mat + transpose(conjg(mat))

    allocate(egnvecs(n,n))
    allocate(iwork(nmax))
    allocate(egnvls(n), rwork(nmax))
    allocate(work(nmax))
    allocate(egndiag(n,n),mat2(n,n),hermmat1(n,n),Q(n,n),diff(n,n),matcheck(n,n))

    if(count((hermmat==transpose(conjg(hermmat))))==size(hermmat)) then
        print *, 'Generated matrix is Hermitian'
    else
        print *, 'Error - Generated matrix is not Hermitian'
        stop
    end if

    egnvecs = hermmat !'zheevd' function reduces input matrix to eigenvector matrix
    !query workspace
    lwork = -1
    liwork = -1
    lrwork = -1
    call zheevd('v','u',n,egnvecs,n,egnvls,work,lwork,rwork,lrwork,iwork,liwork,info)
    ! lwork = min(nmax,int(work(1)))
    ! lrwork = min(nmax,int(rwork(1)))
    ! liwork = min(nmax,iwork(1))
    lwork = int(work(1))
    lrwork = int(rwork(1))
    liwork = iwork(1)
    !eigenvalues & eigenvectors
    call zheevd('v','u',n,egnvecs,n,egnvls,work,lwork,rwork,lrwork,iwork,liwork,info)
    !check
    if (info>0) then
        print *, 'Error from LAPACK subroutine'
        stop
    end if

    !checking if eigenvalues & eigenvectors are correct or not
    Q = transpose(conjg(egnvecs))
    call cbasis_change(hermmat,Q,hermmat1)
    egndiag = 0 !creating matrix with eigenvalues as diagonal elements
    do i=1,size(egnvls)
        egndiag(i,i) = egnvls(i)
    end do
    diff = hermmat1-egndiag
    if(count(abs(diff)<error)==size(diff)) then
        print *, 'Correct eigenvalues & vectors'
    else
        print *, 'Error - Incorrect eigenvalues & vectors'
        stop
    end if

    !storing eigenvalues & eigenvectors on binary files
    open(newunit=io1, file="eigenvalues1.bin", status="replace", action="write", position="rewind", iostat=info)
    if(info/=0) stop "Cannot open file-error"
    write(io1,*) egnvls
    close(io1)
    open(newunit=io2, file="eigenvectors1.bin", status="replace", action="write", position="rewind", iostat=info)
    if(info/=0) stop "Cannot open file-error"
    do i=1,size(egnvecs,1)
        write(io2,*) egnvecs(i,:)
    end do
    close(io2)

    !generating another random matrix
    call random_mtx(n,mat1)
    mat1 = mat1*100
    if((count(mat1==mat))==size(mat)) then
        print *, 'Error- 2 same random matrices generated'
        stop
    else
        print *, 'No error- 2 different random matrices generated'
    end if

    !changing basis of 2nd random matrix with eigenvectors of 1st hermitian matrix
    call cbasis_change(mat1,Q,mat2)
    !checking whether basis change is correct or not
    Q = transpose(conjg(Q))
    call cbasis_change(mat2,Q,matcheck)
    diff = matcheck-mat1
    if(count(abs(diff)<error)==size(diff)) then
        print *, 'Same matrix obtained after reversing basis change - No error'
    else
        print *, 'Error - Different matrix obtained after reversing basis change'
    end if

end program main1