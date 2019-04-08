module integral_mod
    implicit none
    real(8)::i1=2.15888308336d0,i2=0.74682413281234d0,i3=0.7266423406817255d0
contains
    function f1(x)
        real(8) :: x,f1
        f1 = dlog(x)
    end function f1
    function f2(x)
        real(8) :: x,f2
        f2 = dexp(-x*x)
    end function f2
    function f3(x)
        real(8) :: x,f3
        f3 = 1.d0/(x*x+1.d0)
    end function f3
    
    function Trapzoied(x0,xn,n,f)
        real(8)::x,x0,xn,h,f,Trapzoied
        integer(16)::n
        h = (xn-x0)/n
        x = x0 + h
        Trapzoied = (f(x0)+f(xn))/2
        do while(x<xn)
            Trapzoied = Trapzoied + f(x)
            x = x + h
        end do
        Trapzoied = Trapzoied * h
    end function Trapzoied

    function Simpson(x0,xn,n,f)
        real(8)::x,x0,xn,h,f,Simpson, Sodd, Seven
        integer(16)::n
        h = (xn-x0)/n
        x = x0 + 2*h
        Simpson = f(x0)+f(xn)
        Sodd = f(x-h)
        Seven = 0.d0
        do while(x<xn)
            Sodd = Sodd + f(x+h)
            Seven = Seven + f(x)
            x = x + 2*h
        end do
        Simpson = (Simpson + Sodd*4.d0 + Seven*2.d0)*h/3.d0
    end function Simpson
end module integral_mod

program main
    use integral_mod
    ! implicit none
    real(8)::T,S,x0,xn,err,err1
    integer(16)::n
    
    x0 = 2.0d0
    xn = 4.0d0
    n=100
    err=1
    do while (err>0.0000001d0)
        T = Trapzoied(x0,xn,n,f1)
        err = abs(T-i1)
        n=n+1
    enddo
    write(*,*)'Trapzoied: ', n-1,err

    x0 = 2.0d0
    xn = 4.0d0
    n=10
    err=1
    do while (err>0.0000001d0)
        S = Simpson(x0,xn,n,f1)
        err = abs(S-i1)
        n=n+2
    enddo
    write(*,*)'Simpson: ', n-2,err

    x0 = 0.0d0
    xn = 1.0d0
    n=100
    err=1
    do while (err>0.000000001d0)
        T = Trapzoied(x0,xn,n,f2)
        err = abs(T-i2)
        n=n+1
    enddo
    write(*,*)'Trapzoied: ', n-1,err

    x0 = 0.d0
    xn = 1.d0
    n=10
    err=1
    do while (err>0.000000001d0)
        T = Simpson(x0,xn,n,f2)
        err = abs(T-i2)
        n=n+2
    enddo
    write(*,*)'Simpson: ', n-2,err

    x0 = 0.5d0
    xn = 2.5d0
    n=100
    err=1
    do while (err>0.00000000001d0)
        T = Trapzoied(x0,xn,n,f3)
        err = abs(T-i3)
        n=n+1
    enddo
    write(*,*)'Trapzoied: ', n-1,err

    x0 = 0.5d0
    xn = 2.5d0
    n=10
    err=1
    do while (err>0.00000000001d0)
        T = Simpson(x0,xn,n,f3)
        err = abs(T-i3)
        n=n+2
    enddo
    write(*,*)'Simpson: ', n-2,err

end program main





