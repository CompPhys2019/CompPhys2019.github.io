program main
    implicit none
    real(8),allocatable:: a(:,:),eigvals(:),r
    real,allocatable:: a1(:,:),eigvals1(:),r1

    integer::i,j,n

    open(10,file='ratio.dat')
    do n = 2, 101
        allocate(a(n,n),a1(n,n),eigvals(n),eigvals1(n))

        do i = 1, n
            do j = 1, n
                a(i,j) = 1.d0/(i+j-1)
                a1(i,j) = 1.0/(i+j-1)
            end do
        end do

        call diag(n,a,eigvals)
        eigvals = dabs(eigvals)
        r = dlog(maxval(eigvals)/minval(eigvals))
    
        call sdiag(n,a1,eigvals1)
        eigvals1 = abs(eigvals1)
        r1 = log(maxval(eigvals1)/minval(eigvals1))

        ! write(*,*)eigvals
        deallocate(a,a1,eigvals,eigvals1)
        write(10,*)n,r,r1
    end do

end program main

subroutine diag(nd,vec,val)
    implicit none
    integer::nd,lwork,info
    real(8)::vec(nd,nd),val(nd),work(nd*(3+nd/2))

    lwork = nd*(3+nd/2)

    call dsyev('V','U',nd,vec,nd,val,work,lwork,info)
    
end subroutine diag

subroutine sdiag(nd,vec,val)
    implicit none
    integer::nd,lwork,info
    real::vec(nd,nd),val(nd),work(3*nd)

    lwork = 3*nd

    call ssyev('V','U',nd,vec,nd,val,work,lwork,info)
    
end subroutine sdiag
