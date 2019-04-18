module det_mod
    implicit none
contains

recursive function det(matrix) result(laplace_det)
    real(8) :: matrix(:,:)
    integer :: i, n, p, q
    real(8) :: laplace_det, det1
    real(8), allocatable :: cf(:,:)

    n = size(matrix, 1) 

    if (n == 1) then  
        det1 = matrix(1,1)
    else
        det1 = 0
        do i = 1, n  
            allocate( cf(n-1, n-1) )
            cf = cofactor( matrix, i, 1 )

            det1 = det1 + ((-1)**(i+1))* matrix( i, 1 ) * det( cf )
            deallocate(cf)
        end do
    end if

    laplace_det = det1

end function

function cofactor(matrix, mI, mJ)
  real(8), dimension(:,:) :: matrix
  integer :: mI, mJ
  integer :: msize(2), i, j, k, l, n
  real(8), dimension(:,:), allocatable :: cofactor
  msize = shape(matrix)
  n = msize(1)

  allocate(cofactor(n-1, n-1))
  l=0
  k = 1
  do i=1, n
   if (i .ne. mI) then
     l = 1
     do j=1, n
       if (j .ne. mJ) then
         cofactor(k,l) = matrix(i,j)
         l = l+ 1
       end if
     end do
     k = k+ 1
   end if
  end do
return
end function

end module

program main
    use det_mod
    implicit none
    integer::i
    real(8):: a(3,3),b(3,3),c(3,3),e(3,3)
    real(8):: d(3,3),d1(3,3),d2(3,3),d3(3,3),d4(3,3),d5(3,3),d6(3,3),d7(3,3),d8(3,3)
    real(8):: det1, det2

    a(:,1) = [1.d0,1.d0,-1.d0]
    a(:,2) = [1.d0,2.d0,-2.d0]
    a(:,3) = [-2.d0,1.d0,1.d0]

    b(:,1) = [0.d0,1.d0,2.d0]
    b(:,2) = [1.d0,2.d0,1.d0]
    b(:,3) = [2.d0,1.d0,0.d0]

    c(:,1) = [3.d0,2.d0,1.d0]
    c(:,2) = [1.d0,2.d0,3.d0]
    c(:,3) = [1.d0,-1.d0,1.d0]

    write(*,*)"A = "
    write(*,'(3f8.2)')a
    write(*,*)' '
    write(*,*)"B = "
    write(*,'(3f8.2)')b
    write(*,*)' '
    write(*,*)"C = "
    write(*,'(3f8.2)')c

    write(*,*)"--------------------------"

    d = matmul(matmul(a,b),c)
    d1 = matmul(a,matmul(b,c))

    write(*,*)"ABC = "
    write(*,'(3f8.2)')d
    write(*,*)' '
    write(*,*)"A(BC) = "
    write(*,'(3f8.2)')d1

    write(*,*)"--------------------------"

    d2 = matmul(a,b+c)
    d3 = matmul(a,b) + matmul(a,c)

    write(*,*)"A(B+C) = "
    write(*,'(3f8.2)')d2
    write(*,*)' '
    write(*,*)"A(BC) = "
    write(*,'(3f8.2)')d3

    write(*,*)"--------------------------"

    d4 = transpose(matmul(a,b))
    d5 = matmul(transpose(b),transpose(a))

    write(*,*)"（AB)^T = "
    write(*,'(3f8.2)')d4
    write(*,*)' '
    write(*,*)"B^T A^T = "
    write(*,'(3f8.2)')d5

    write(*,*)"--------------------------"

    call matrixinv(matmul(a,b),d6,3)
    call matrixinv(a,d7,3)
    call matrixinv(b,d8,3)

    write(*,*)"（AB)^-1 = "
    write(*,'(3f8.2)')d6
    write(*,*)' '
    write(*,*)"B^-1 A^-1 = "
    write(*,'(3f8.2)')matmul(d8,d7)

    a(:,1) = [1.d0,0.d0, 0.d0]
    a(:,2) = [0.d0,2.d0, 0.d0]
    a(:,3) = [0.d0,0.d0,2.d0]

    det1 = det(matmul(a,b))
    det2 = det(a)*det(b)
    write(*,*)'det(a,b) = ',det1
    write(*,*)''
    write(*,*)'det(a)det(b) =', det2
end program main

subroutine matrixinv(a,b,n)
    ! subroutine to calculate the inverse of a matrix using Gauss-Jordan elimination
    ! the inverse of matrix a(n,n) is calculated and stored in the matrix b(n,n)
    ! input ...
    ! a(n,n) - array of coefficients for matrix A
    ! n      - dimension
    ! output ...
    ! c(n,n) - inverse matrix of A
    integer :: i,j,k,l,m,n,irow
    real(8):: big,a(n,n),b(n,n),dum
    !build the identity matrix
    do i = 1,n
    do j = 1,n
    b(i,j) = 0.0
    end do
    b(i,i) = 1.0
    end do
    do i = 1,n ! this is the big loop over all the columns of a(n,n)
    ! in case the entry a(i,i) is zero, we need to find a good pivot; this pivot
    ! is chosen as the largest value on the column i from a(j,i) with j = 1,n
    big = a(i,i)
    do j = i,n
    if (a(j,i).gt.big) then
    big = a(j,i)
    irow = j
    end if
    end do
    ! interchange lines i with irow for both a() and b() matrices
    if (big.gt.a(i,i)) then
    do k = 1,n
    dum = a(i,k) ! matrix a()
    a(i,k) = a(irow,k)
    a(irow,k) = dum
    dum = b(i,k) ! matrix b()
    b(i,k) = b(irow,k)
    b(irow,k) = dum
    end do
    end if
    ! divide all entries in line i from a(i,j) by the value a(i,i);
    ! same operation for the identity matrix
    dum = a(i,i)
    do j = 1,n
    a(i,j) = a(i,j)/dum
    b(i,j) = b(i,j)/dum
    end do
    ! make zero all entries in the column a(j,i); same operation for indent()
    do j = i+1,n
    dum = a(j,i)
    do k = 1,n
    a(j,k) = a(j,k) - dum*a(i,k)
    b(j,k) = b(j,k) - dum*b(i,k)
    end do
    end do
    end do
    ! substract appropiate multiple of row j from row j-1
    do i = 1,n-1
    do j = i+1,n
    dum = a(i,j)
    do l = 1,n
    a(i,l) = a(i,l)-dum*a(j,l)
    b(i,l) = b(i,l)-dum*b(j,l)
    end do
    end do
    end do
end subroutine matrixinv
