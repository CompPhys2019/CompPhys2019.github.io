module random
    ! Description =================================================================================
    !   rn is an overloaded function, which rn() generates real random numbers from 0 to 1,
    !   not including 1, and rn(n1,n2) generates integer numbers from n1 to n2.
    ! ---------------------------------------------------------------------------------------------

    ! Arguments -----------------------------------------------------------------------------------
    !   seed
    !     integer type, initial seed for initiate random generator
    !   a_min, a_max
    !     integer type, minimal and maximal random number
    ! ---------------------------------------------------------------------------------------------

    ! Function rn() ------------------------------------------------------------------------------
    !   generates real random numbers from 0 to 1, not including 1
    ! ---------------------------------------------------------------------------------------------

    ! Function rn(a_min,a_max) -------------------------------------------------------------------
    !   generates integer random numbers from a_min to a_max, including a_min and a_max
    ! ---------------------------------------------------------------------------------------------

    ! Subroutine initrn(seed) --------------------------------------------------------------------
    !   initate random generator
    ! ---------------------------------------------------------------------------------------------

    ! =============================================================================================

    interface rn
        module procedure ran_real,ran_int
    end interface
contains
    subroutine initrn(seed)
        implicit none
        integer(8) :: irmax
        real(8)    :: dmu64
        integer(8) :: ran64,mul64,add64
        common/bran64/dmu64,ran64,mul64,add64

        integer::seed

        irmax=2_8**31
        irmax=2*(irmax**2-1)+1
        mul64=2862933555777941757_8
        add64=1013904243
        dmu64=0.5d0/dble(irmax)

        ran64=abs((seed*mul64)/5+5265361)
    end subroutine initrn

    function ran_real()
        implicit none
        real(8)::ran_real
        real(8)    :: dmu64
        integer(8) :: ran64,mul64,add64
        common/bran64/dmu64,ran64,mul64,add64

        ran64=ran64*mul64+add64
        ran_real=0.5d0+dmu64*dble(ran64)
    end function ran_real

    function ran_int(a_min,a_max)
        integer::a_min,a_max,ran_int
        ran_int=int((a_max-a_min+1)*ran_real())+a_min
    end function ran_int
end module random