program main
    implicit none
    real(8)::a(24,5)
    integer::i
    open(10,file='temp.dat')

    do i = 1, 24
        read(10,*)a(i,:)
    end do 

    do i = 1, 24
        write(*,'(5f8.4)')a(i,:)
    end do 
end program main