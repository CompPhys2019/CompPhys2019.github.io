! #################################################################################################
! Random Walk Algorithm
! AUTHOR: W.J.Zhu
! Date: May 10, 2019
! Mail: wjzhu@mail.bnu.edu.cn
! #################################################################################################

! *************************************************************************************************
! This is Fortran90 implementation of Random Walks, including simple RW, non-reversal RW and 
! self-avoiding RW.

! Algorithm: 
!   For simple RW, there is no restriction to each random step.
!   For non-reversal RW, backword step is not allowed.
!   For self-avoiding RW, no lattice size can be visited more than onece

! Measurement:
!   <x(N)>, <y(N)>, <\Delta R^2(N)> = <x^2(N)> - <x(N)>^2 + <y^2(N)> -<y(N)>^2
! *************************************************************************************************

module random
    ! Description =================================================================================
    !   rn is an random function, which rn() generates real random numbers from 0 to 1,
    !   not including 1, 
    !   rnint generates integer numbers from n1 to n2.
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

    ! Function rnint(a_min,a_max) -------------------------------------------------------------------
    !   generates integer random numbers from a_min to a_max, including a_min and a_max
    ! ---------------------------------------------------------------------------------------------

    ! Subroutine initrn(seed) --------------------------------------------------------------------
    !   initate random generator
    ! ---------------------------------------------------------------------------------------------

    ! =============================================================================================

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

    function rn()
        implicit none
        real(8)::rn
        real(8)    :: dmu64
        integer(8) :: ran64,mul64,add64
        common/bran64/dmu64,ran64,mul64,add64

        ran64=ran64*mul64+add64
        rn=0.5d0+dmu64*dble(ran64)
    end function rn

    function rnint(a_min,a_max)
        integer::a_min,a_max,rnint
        rnint=int((a_max-a_min+1)*rn())+a_min
    end function rnint
end module random

module RW
    use random
    implicit none
    integer::move(4,2)
contains
    subroutine init_RW
        implicit none
        ! For each step, there are four possible moving steps, moving left, right, down or up.
        ! 1: (x,y) --> (x-1,y)
        ! 2: (x,y) --> (x+1,y)
        ! 3: (x,y) --> (x,y-1)
        ! 4: (x,y) --> (x,y+1)
        move(:,1) = [-1,1,0,0]
        move(:,2) = [0,0,-1,1]
    end subroutine init_RW

    subroutine Simple_RW(r)
        ! Simple RW Algorithm only needs current position to update. 
        implicit none
        integer::r(2),i

        i = rnint(1,4) ! choose move step
        r = r + move(i,:)
    end subroutine Simple_RW

    subroutine Non_Reversal_RW(r,r_pre)
        ! Previous position needed
        implicit none
        integer::r(2),r_pre(2),r1(2)
        integer::i,j,k

        r1 = r
        ! ruolette wheel selction
        ! Suppose (x+1,y) is previous position, we should select 1,3,4 from index i
        ! Sequencely, we can assign 1-->1, 2-->3, 3-->4. 
        i = rnint(1,3) ! choose move step
        k = 1
        do j = 1, 4
            if(r(1)+move(j,1)==r_pre(1) .and. r(2)+move(j,2)==r_pre(2))cycle
            if(k==i)then
                r = r + move(j,:)
                exit
            end if
            k = k + 1
        end do

        r_pre = r1
    end subroutine Non_Reversal_RW

    subroutine Self_Avoiding_RW(r,mask,status)
        ! mask is an array labeling visited sites, false for unvisited and true for visited
        use random
        implicit none
        integer::r(2)
        logical::status
        logical,allocatable::mask(:,:)
        integer::i,j,cnt,neib(4)

        status = .false.
        mask(r(1),r(2))= .true.

        ! count unvisited neighbor sites
        cnt = 0
        do i = 1, 4 
            if(mask(r(1)+move(i,1),r(2)+move(i,2)))cycle
            cnt = cnt + 1 
            neib(cnt) = i
        end do

        if (cnt==0) then
            status = .true.
        else
            i = rnint(1,cnt)
            j = neib(i)
            r = r + move(j,:)
        endif

    end subroutine Self_Avoiding_RW
end module RW

module measure_mod
    implicit none
    integer::Nsam, Nbin
    real(8)::rm(2),dr2(2)
contains

    subroutine init_measure
        implicit none
        rm = 0.d0
        dr2 = 0.d0
        
    end subroutine init_measure

    subroutine measure(r)
        implicit none
        integer::r(2)

        rm = rm + r
        dr2 = dr2 + r**2

    end subroutine measure

    subroutine results
        implicit none
        
        rm = rm / Nsam
        dr2 = dr2 / Nsam
        dr2 = dr2 - rm**2 

    end subroutine results
end module measure_mod

program main
    use random
    use measure_mod
    use RW
    implicit none
    integer::N, r(2), r_pre(2), seed
    ! N: walk steps
    ! x,y: coordinate

    logical,allocatable::mask(:,:) ! labels visited sites as true. 
    logical::status

    integer::i,j,k

    N = 6 ! set default
    allocate(mask(-N:N,-N:N))

    read(*,*)Nsam, Nbin, seed

    write(*,*)'Nsam = ', Nsam, ',    Nbin = ', Nbin
    write(*,*)'------------------------------------------------------------------------'

    call initrn(seed)

    call init_RW

    write(*,*)'Simple RW: '
    open(10,file='SRW.dat',position='append',status='unknown')
    do k = 1, Nbin
        call init_measure
        do i = 1, Nsam
            r = 0
            do j = 1, N
                call Simple_RW(r)
            end do
            call measure(r)
        end do
        call results
        write(*,*)rm,sum(dr2)
        write(10,*)rm,sum(dr2)
    end do
    close(10)

    write(*,*)'Non Reversal RW: '
    open(10,file='NRRW.dat',position='append',status='unknown')
    do k = 1, Nbin
        call init_measure
        do i = 1, Nsam
            r = 0
            call Simple_RW(r)
            r_pre = 0
            do j = 1, N-1
                call Non_Reversal_RW(r,r_pre)
            end do
            call measure(r)
        end do
        call results
        write(*,*)rm,sum(dr2)
        write(10,*)rm,sum(dr2)
    end do
    close(10)

    write(*,*)'Self Avoiding RW: '
    do k = 1, Nbin
        call init_measure
        do i = 1, Nsam
            r = 0
            mask = .false.
            do j = 1, N
                call Self_Avoiding_RW(r,mask,status)
                if(status) exit
            end do
            call measure(r)
        end do
        call results
        write(*,*)rm,sum(dr2)
        write(10,*)rm,sum(dr2)
    end do
    close(10)

end program main



