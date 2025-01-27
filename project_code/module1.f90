module mymodule1 !MODULE HAVING SUBROUTINES TO PRINT MATRIX, GENERATE RANDOM COMPLEX MATRIX & CHANGE BASIS
    use, intrinsic :: iso_fortran_env, only: a=>real64
    implicit none
    
    private
    public rprint_mat, cprint_mat, random_mtx, rbasis_change, cbasis_change

contains
    subroutine rprint_mat(B)
        real(a), intent(in) :: B(:,:)
        integer :: i

        do i = 1, size(B,1)
            print *, B(i,:)            
        end do
    end subroutine rprint_mat

    subroutine cprint_mat(B)
        complex(a), intent(in) :: B(:,:)
        integer :: i

        do i = 1, size(B,1)
            print *, B(i,:)            
        end do
    end subroutine cprint_mat

    subroutine random_mtx(n,C)
        integer, intent(in) :: n
        integer :: i, j, clock
        integer, allocatable :: seed(:)
        real(a), allocatable :: realmat(:,:), imgmat(:,:)
        complex(a), allocatable, intent(out) :: C(:,:)

        allocate(realmat(n,n),imgmat(n,n))
        allocate(C(n,n))

        call random_seed(size=j)
        allocate(seed(j))
        call system_clock(count=clock)
        seed = clock + 37 * (/(i-1, i=1, j)/)
        call random_seed(put=seed)
        deallocate(seed)

        call random_number(realmat)
        call random_number(imgmat)

        loop1: do i=1,size(C,1)
            loop2: do j=1,size(C,2)
                C(i,j) = complex(realmat(i,j),imgmat(i,j))
            end do loop2
        end do loop1

    end subroutine random_mtx

    subroutine rbasis_change(input,bc_mat,output)
        real(a),intent(in) :: input(:,:), bc_mat(:,:)
        real(a),intent(out) :: output(:,:)
        real(a), allocatable :: temp(:,:)

        allocate(temp(size(input,1),size(input,2)))

        temp = matmul(input,transpose(bc_mat))
        output = matmul(bc_mat,temp)

    end subroutine rbasis_change

    subroutine cbasis_change(input,bc_mat,output)
        complex(a),intent(in) :: input(:,:), bc_mat(:,:)
        complex(a),intent(out) :: output(:,:)
        complex(a), allocatable :: temp(:,:)

        allocate(temp(size(input,1),size(input,2)))

        temp = matmul(input,transpose(conjg(bc_mat)))
        output = matmul(bc_mat,temp)

    end subroutine cbasis_change
end module mymodule1