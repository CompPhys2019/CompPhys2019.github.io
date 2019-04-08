program main
    implicit none
    integer::i
    integer:: a(3,3),b(3,3),c(3,3)
    integer:: d(3,3),d1(3,3),d2(3,3),d3(3,3),d4(3,3),d5(3,3),d6(3,3)

    a(:,1) = [1,1,-1]
    a(:,2) = [1,2,-2]
    a(:,3) = [-2,1,1]

    b(:,1) = [0,1,2]
    b(:,2) = [1,2,1]
    b(:,3) = [2,1,0]

    c(:,1) = [3,2,1]
    c(:,2) = [1,2,3]
    c(:,3) = [1,-1,1]

    write(*,*)"A = "
    write(*,'(3I4)')a
    write(*,*)' '
    write(*,*)"B = "
    write(*,'(3I4)')b
    write(*,*)' '
    write(*,*)"C = "
    write(*,'(3I4)')c

    write(*,*)"--------------------------"

    d = matmul(matmul(a,b),c)
    d1 = matmul(a,matmul(b,c))

    write(*,*)"ABC = "
    write(*,'(3I4)')d
    write(*,*)' '
    write(*,*)"A(BC) = "
    write(*,'(3I4)')d1

    write(*,*)"--------------------------"

    d2 = matmul(a,b+c)
    d3 = matmul(a,b) + matmul(a,c)

    write(*,*)"A(B+C) = "
    write(*,'(3I4)')d2
    write(*,*)' '
    write(*,*)"A(BC) = "
    write(*,'(3I4)')d3

    write(*,*)"--------------------------"

    d4 = transpose(matmul(a,b))
    d5 = matmul(transpose(b),transpose(a))

    write(*,*)"（AB)^T = "
    write(*,'(3I4)')d4
    write(*,*)' '
    write(*,*)"A^T B^T = "
    write(*,'(3I4)')d5

    write(*,*)"--------------------------"

    write(*,'(3f8.4)')zInverse(3,real(a,8))

contains
    function zInverse(n, a)  result(ra)

    integer::n,lda,ipiv(n),info,lwork

    real(8)::a(n,n),ra(n,n),work(n)

      ra=a

      lwork=n

      lda=n

      call zgetrf(n, n, ra, lda, ipiv, info)

      if(info/=0) write(0,*) 'Error occured in zgetrf!'

      call zgetri(n, ra, lda, ipiv, work, lwork, info)

      if(info/=0) write(0,*) 'Error occured in zgetri!'

    end function

end program main
