module diff
contains
    function dfun(h,c,n)
        implicit none
        real(8),allocatable::h(:),dfun(:)
        real(8)::c(-3:3)
        integer::i,n

        dfun = 0.d0*h
        do i = -3, 3
            dfun = dfun + c(i)*exp(i*h)/(sin(i*h)**3+cos(i*h)**3)
        end do
        dfun = dfun/h**n
    end function dfun
end module diff

program main
    use diff
    implicit none
    real(8)::c(-3:3)
    real(8),allocatable::h(:),df(:)
    integer::n,i

    n = 10
    allocate(h(n),df(n))
    h = [(1.d0/2**i,i=1,n)]

    ! two-point forward first
    c = 0.d0
    c(0:1) = [-1,1]
    df = dfun(h,c,1)
    write(*,*)dfun(h,c,1)
    write(*,*)''


    ! three-point first
    c = 0
    c(-1:1)=[-0.5,0.0,0.5] 
    write(*,*)dfun(h,c,1)
    write(*,*)''

    ! ! five-point first
    c = 0
    c(-2:2)=[1.0/12,-2.0/3,0.0,2.0/3,-1.0/12] 
    write(*,*)dfun(h,c,1)
    write(*,*)''

    ! ! second derivative
    c = 0
    c(-1:1)=[1.0,-2.0,1.0]
    write(*,*)dfun(h,c,2)
    write(*,*)''


end program main

subroutine inverse(a,c,n)
!============================================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! the original matrix a(n,n) will be destroyed 
! during the calculation
!===========================================================
implicit none 
integer n
double precision a(n,n), c(n,n)
double precision L(n,n), U(n,n), b(n), d(n), x(n)
double precision coeff
integer i, j, k

! step 0: initialization for matrices L and U and b
! Fortran 90/95 aloows such operations on matrices
L=0.0
U=0.0
b=0.0

! step 1: forward elimination
do k=1, n-1
   do i=k+1,n
      coeff=a(i,k)/a(k,k)
      L(i,k) = coeff
      do j=k+1,n
         a(i,j) = a(i,j)-coeff*a(k,j)
      end do
   end do
end do

! Step 2: prepare L and U matrices 
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
do i=1,n
  L(i,i) = 1.0
end do
! U matrix is the upper triangular part of A
do j=1,n
  do i=1,j
    U(i,j) = a(i,j)
  end do
end do

! Step 3: compute columns of the inverse matrix C
do k=1,n
  b(k)=1.0
  d(1) = b(1)
! Step 3a: Solve Ld=b using the forward substitution
  do i=2,n
    d(i)=b(i)
    do j=1,i-1
      d(i) = d(i) - L(i,j)*d(j)
    end do
  end do
! Step 3b: Solve Ux=d using the back substitution
  x(n)=d(n)/U(n,n)
  do i = n-1,1,-1
    x(i) = d(i)
    do j=n,i+1,-1
      x(i)=x(i)-U(i,j)*x(j)
    end do
    x(i) = x(i)/u(i,i)
  end do
! Step 3c: fill the solutions x(n) into column k of C
  do i=1,n
    c(i,k) = x(i)
  end do
  b(k)=0.0
end do
end subroutine inverse