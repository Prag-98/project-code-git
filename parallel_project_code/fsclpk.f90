program fsclpk
    use, intrinsic :: iso_fortran_env, only: a=>real64, b=>int64
    implicit none

    complex(a), allocatable :: hermmat(:,:), check(:,:), diff(:,:)
    complex(a), allocatable :: egnvecs(:,:), work(:),work1(:)
    real(a), allocatable :: egnvls(:), rwork(:)
    integer(b), allocatable :: iwork(:)
    complex(a), allocatable :: Q(:,:), temp(:,:), diag(:,:), hermmat1(:,:)
    complex(a), allocatable :: mat2(:,:), mat2nb(:,:), writemat(:,:)
    
    integer(b) :: desca(9), descb(9)

    integer :: n, nb, nprows, npcols, nbmx
    integer :: cntxt, myrow, mycol, info, ierr
    integer :: lwork, liwork, lrwork, info1
    integer :: i, j, io, io1, io2

    real(a) :: compare, error

    complex(a) :: alpha, beta

    !n = 10
    !n = 1000
    !n = 10000
    !n = 29500
    n = 13056
    
    !nb = 2
    nb = 32

    nprows=2
    npcols=2

    !nbmx = 6
    !nbmx = 512
    !nbmx = 5008
    !nbmx = 14752
    nbmx = 6528
    
    error = 1e-6*1.0_a

    allocate(hermmat(nbmx,nbmx), check(nbmx,nbmx), diff(nbmx,nbmx))
    allocate(egnvecs(nbmx,nbmx), egnvls(n))
    allocate(Q(nbmx,nbmx), temp(nbmx,nbmx), diag(nbmx,nbmx), hermmat1(nbmx,nbmx))
    allocate(work1(nb))
    allocate(mat2(nbmx,nbmx), mat2nb(nbmx,nbmx))

    check = dcmplx(0.0,0.0)
    alpha = dcmplx(1.0,0.0)
    beta = dcmplx(1.0,0.0)
    Q = dcmplx(0.0,0.0)
    temp = dcmplx(0.0,0.0)
    hermmat1 = dcmplx(0.0,0.0)
    diag = dcmplx(0.0,0.0)

    !intializing BLACS
    call blacs_get(0, 0, cntxt)
    !intializing process grid
    call blacs_gridinit(cntxt, 'R', nprows, npcols)
        call blacs_gridinfo(cntxt, nprows, npcols, myrow, mycol)

        !array descriptor for input hermitian matrix & eigenvectors matrix
        call descinit(desca, n, n, nb, nb, 0, 0, cntxt, nbmx, info)
        if ( info/=0 ) then
            if ( myrow==0.AND.mycol==0 ) then
                print *, 'Error in assigning array descriptor', info
                call blacs_abort(cntxt, ierr)
            end if
        end if
        call descinit(descb, n, n, n, n, 0, 0, cntxt, n, info)
        if ( info/=0 ) then
            print *, 'Error in assigning array descriptor'
            call blacs_abort(cntxt, ierr)
        end if
        if ( myrow==0.AND.mycol==0 ) then
            allocate(writemat(n,n))
        end if
        
        !creating a hermitian matrix
        loop1: do i = 1, n
            loop2: do j = 1, n
                if ( i==j ) then
                    call pzelset(hermmat, i, j, desca, dcmplx((1.0_a/((i+j-1)*1.0_a))+1.0_a,0.0))
                else
                    call pzelset(hermmat, i, j, desca, dcmplx((1.0_a/((i+j-1)*1.0_a)),(i-j)))
                end if
            end do loop2
        end do loop1

        !checking if matrix is hermitian
        call pztranc(n, n, alpha, hermmat, 1, 1, desca, beta, check, 1, 1, desca)
        diff = check-hermmat
        call zgsum2d(cntxt, 'All', 'I-ring', nbmx, nbmx, diff, nbmx, 0, 0)
        if ( myrow==0.AND.mycol==0 ) then
            compare = zabs(sum(diff))
            if ( compare>error ) then
                print *, '1st Matrix is not Hermitian', compare
                call blacs_abort(cntxt, ierr)
            else
                print *, '1st Matrix is Hermitian'
            end if
        end if

        check = hermmat

        !query workspace
        lwork=-1
        lrwork=-1
        liwork=-1
        allocate(work(1),rwork(1),iwork(1))
        call pzheevd('V', 'U', n, hermmat, 1, 1, desca, egnvls, egnvecs, 1, 1, &
                     desca, work, lwork, rwork, lrwork, iwork, liwork, info1)
        if ( info1/=0 ) then
            if ( myrow==0.AND.mycol==0 ) then
                print *, 'Error from PZHEEVD routine', info1
                call blacs_abort(cntxt, ierr)
            end if
        end if
        lwork = int(work(1))
        liwork = int(iwork(1))
        lrwork = int(rwork(1))
        deallocate(work,rwork,iwork)
        allocate(work(lwork),rwork(lrwork),iwork(liwork))
        !find eigenvalues & eigenvectors
        call pzheevd('V', 'U', n, hermmat, 1, 1, desca, egnvls, egnvecs, 1, 1, &
                     desca, work, lwork, rwork, lrwork, iwork, liwork, info1)
        if ( info1/=0 ) then
            if ( myrow==0.AND.mycol==0 ) then
                print *, 'Error from PZHEEVD routine', info1
                call blacs_abort(cntxt, ierr)
            end if
        end if

        !checking if eigenvalues & vectors are correct or not
        loop3: do i = 1, n
            loop4: do j = 1, n
                if ( i==j ) then
                    call pzelset(diag, i, j, desca, dcmplx(egnvls(i),0.0))
                else
                    call pzelset(diag, i, j, desca, dcmplx(0.0,0.0))
                end if
            end do loop4
        end do loop3
        call pztranc(n, n, alpha, egnvecs, 1, 1, desca, beta, Q, 1, 1, desca)
        call pzgemm('N', 'C', n, n, n, alpha, check, 1, 1, desca, Q, 1, 1, desca, &
                    beta, temp, 1, 1, desca)
        call pzgemm('N', 'N', n, n, n, alpha, Q, 1, 1, desca, temp, 1, 1, desca, &
                    beta, hermmat1, 1, 1, desca)
        diff = hermmat1-diag
        call zgsum2d(cntxt, 'All', 'I-ring', nbmx, nbmx, diff, nbmx, 0, 0)
        if ( myrow==0.AND.mycol==0 ) then
            compare = zabs(sum(diff))
            if ( compare>error ) then
                print *, 'Eigenvalues & vectors are incorrect', compare
                call blacs_abort(cntxt, ierr)
            else
                print *, 'Eigenvalues & vectors are correct'
            end if
        end if

        !Creating another matrix
        loop5: do i = 1, n
            loop6: do j = 1, n
                call pzelset(mat2, i, j, desca, dcmplx((i*j*1.0_a),((i*1.0_a)/j)))
            end do loop6
        end do loop5

        !checking if matrix is hermitian
        check = dcmplx(0.0,0.0)
        call pztranc(n, n, alpha, mat2, 1, 1, desca, beta, check, 1, 1, desca)
        diff = check-mat2
        call zgsum2d(cntxt, 'All', 'I-ring', nbmx, nbmx, diff, nbmx, 0, 0)
        if ( myrow==0.AND.mycol==0 ) then
            compare = zabs(sum(diff))
            if ( compare>error ) then
                print *, '2nd Matrix is not Hermitian', compare
            else
                print *, '2nd Matrix is Hermitian'
            end if
        end if

        !changing basis of 2nd matrix
        temp = dcmplx(0.0,0.0)
        mat2nb = dcmplx(0.0,0.0)
        call pzgemm('N', 'C', n, n, n, alpha, mat2, 1, 1, desca, Q, 1, 1, desca, &
                    beta, temp, 1, 1, desca)
        call pzgemm('N', 'N', n, n, n, alpha, Q, 1, 1, desca, temp, 1, 1, desca, &
                    beta, mat2nb, 1, 1, desca)

        !storing eigenvectors & values of 1st matrix and 2nd matrix in new basis in binary files
        call pzgemr2d(n, n, egnvecs, 1, 1, desca, writemat, 1, 1, descb, cntxt)
        if ( myrow==0.AND.mycol==0 ) then
            open(newunit=io, file="eigenvalues2.bin", status="replace", form="unformatted", &
            action="write", position="rewind", iostat=info)
            if ( info/=0 ) then
                print *, 'Error in opening file', info
                call blacs_abort(cntxt, ierr)
            end if
            write(io) egnvls
            close(io)
            open(newunit=io1, file="eigenvectors2.bin", status="replace", form="unformatted", &
            action="write", position="rewind", iostat=info)
            if ( info/=0 ) then
                print *, 'Error in opening file', info
                call blacs_abort(cntxt, ierr)
            end if
            write(io1) writemat
            close(io1)
        end if
        call pzgemr2d(n, n, mat2nb, 1, 1, desca, writemat, 1, 1, descb, cntxt)
        if ( myrow==0.AND.mycol==0 ) then
            open(newunit=io2, file="2ndmatrix.bin", status="replace", form="unformatted", &
            action="write", position="rewind", iostat=info)
            if ( info/=0 ) then
                print *, 'Error in opening file', info
                call blacs_abort(cntxt, ierr)
            end if
            write(io2) writemat
            close(io2)
        end if

    !exiting process grid
    call blacs_gridexit(cntxt)
    !exiting blacs
    call blacs_exit(0)
    
    
end program fsclpk