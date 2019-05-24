! #################################################################################################
! Random Walk Algorithm
! AUTHOR: W.J.Zhu
! Date: May 23, 2019
! Mail: wjzhu@mail.bnu.edu.cn
! #################################################################################################

! =================================================================================================
! This is Fortran90 implementation of simulating molecular dynamics.

! Hints ###########################################################################################
! MD simulation is nothing more than solving 2nd ODE(second order difference equation). The       #
! addtional step is computing acceleration caused by Lennard-Jones potential, which is the most   #
! complicated part of program. Once acceleration calulated, we can follow routines of 2dn ODE.    #
! Contrast to using some familiar algorithm, we use verlet algorithm to simulate 2nd ODE.         # 
!                                                                                                 #
! One should notice that particles locate in box uniformly avoiding divergence.                   #
! #################################################################################################

! Potential ---------------------------------------------------------------------------------------
! The interaction between two particles is the Lennard_Jones potential.
! U(r) = 4 \epsilon [ (\frac \sigma r)^12 - (\frac \sigma r)^6 ]
! f(r) = - \nabla u(r)
! -------------------------------------------------------------------------------------------------

! Periodic ----------------------------------------------------------------------------------------
! The system is specified on square lattice with periodic boundary condition.
! -------------------------------------------------------------------------------------------------

! Algorithm ---------------------------------------------------------------------------------------
! For updating, we choose verlet algorithm
! x_{n+1} = x_n + v_n \Delta t + a_n (\Delta t)^2/2
! v_{n+1} = v_n + (a_{n+1} + a_n ) \Delta t /2
! -------------------------------------------------------------------------------------------------

! Parameters --------------------------------------------------------------------------------------
!   sigma = 3.4d0 * 10^(-10)
!   m = 6.69d0 * 10^(-26)
!   epsilon = 1.65d0 * 10^(-21)
!   k = 1.38d0 * 10^(-23)
! All these are setted to be 1.
! -------------------------------------------------------------------------------------------------

! Measurement -------------------------------------------------------------------------------------
!   E: energy
!   Ek: kinetic energy
!   Ep: potential energy
!   vc: center moment, sum(v)/N
!   T: temperature
!   P: pressure
! -------------------------------------------------------------------------------------------------

! Pseudo-code -------------------------------------------------------------------------------------
! Step 1: initate system
! Step 2: use Verlet algorithm to compute configuration at time t
! Step 3: if t < t_end, continue step 2, otherwise output result
! -------------------------------------------------------------------------------------------------

! =================================================================================================

module random
    ! Description =================================================================================
    !   rn is a random function, which rn() generates real random numbers from 0 to 1,
    ! ---------------------------------------------------------------------------------------------

    ! Arguments -----------------------------------------------------------------------------------
    !   seed
    !     integer type, initial seed for initiate random generator
    ! ---------------------------------------------------------------------------------------------

    ! Function rn() ------------------------------------------------------------------------------
    !   generates real random numbers from 0 to 1, not including 1
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
end module random

module periodic
    implicit none
    integer::L(2)
    ! L(1)--> Lx; L(2)--> Ly. Lx, Ly define system size;

contains
    function disp(r1,r2)
        ! return displacment between r1 and r2. 
        real(8):: r1(2), r2(2), dr(2), disp(2)
        integer::i

        disp = r1 - r2
        do i = 1, 2
            if(disp(i)>0.5d0*L(i))then
                disp(i) = disp(i) - L(i)
            else if (disp(i)<-0.5d0*L(i)) then
                disp(i) = disp(i) + L(i)
            end if
        end do
    end function

    subroutine pbc(r)
        ! re-assign position due to periodic boundary condition
        real(8),allocatable::r(:,:)
        integer:: i, j, N

        N = size(r,dim=1)
        do i = 1, N
            do j = 1, 2
                if (r(i,j)<0.d0) then
                    r(i,j) = r(i,j) + L(j)
                else if (r(i,j)>L(j)) then
                    r(i,j) = r(i,j) - L(j)
                endif
            end do
        end do
    end subroutine

end module periodic

module Lennard_Jones
    use periodic
    implicit none
contains
    function LJ_potential(dr)
        ! Compute LJ potential from displacment between two particle
        real(8)::LJ_potential,dr(2),rm6
        rm6 = 1.d0 / sum(dr*dr)**3
        LJ_potential = 4.d0 * (rm6*rm6-rm6)  
    end function

    function force(dr) ! force(1) --> fx on x direction; force(2) --> fy on y direction.
        ! Compute force as well as acceleration due to m = 1
        ! Note that mass has been setted as 1. So, a = f
        ! f(\vec{r}) = 24 ( 2/r^12 - 1/r^6 ) \vec{r}/r^2
        ! set rm2 = 1/r2, rm6 = rm2^3
        ! f(\vec{r}) = 24 * rm6 * (2*rm6 - 1)*rm2
        implicit none
        real(8)::dr(2), force(2), rm2, rm6, du

        rm2 = 1.d0 / sum(dr*dr)
        rm6 = rm2 * rm2 * rm2

        du = 24*rm6*(2*rm6-1)*rm2
        force = du * dr

    end function force

    subroutine accel(r,a)
        ! compute acceleration given particles' positions r. 
        implicit none
        real(8),allocatable::r(:,:),a(:,:)

        real(8)::dr(2), fr(2)
        integer::i,j, N

        N = size(a,dim=1)
        a = 0.d0        
        do i = 1, N - 1
            do j = i+1, N
                dr = disp(r(i,:),r(j,:))
                fr = force(dr)
                a(i,:) = a(i,:) + fr
                a(j,:) = a(j,:) - fr
            end do
        end do
    end subroutine accel
end module Lennard_Jones

module md
    use random
    use Lennard_Jones
    implicit none
    real(8)::dt,dt2,tend
    ! dt, update interval time; dt2 = dt*dt
    integer::N
    ! L(1)--> Lx; L(2)--> Ly. Lx, Ly define system size; N: particle numbers.

    real(8),allocatable::r(:,:),v(:,:),a(:,:)
    ! allocate(r(N,2),v(N,2),a(N,2))
    ! r: position, r(:,1) --> x; r(:,2) --> y
    ! v: velocity
    ! a: acceleration
    character(32)::fn
contains
    subroutine init_md
        implicit none
        integer::seed
        integer::i,j,k,N1
        real(8)::rx,ry,ratio

        read(*,*)L, N, dt ,tend, seed, fn

        call initrn(seed)

        dt2 = dt * dt

        allocate(r(N,2),v(N,2),a(N,2))

        ratio=sqrt(product(L)/dble(N))

        k = 1
        ry = 0.d0
        do while(ry<L(2))
            rx = 0.d0
            do while(rx<L(1))
                if(k<=N)then
                    r(k,:) = [rx,ry]
                else
                    exit
                endif
                rx = rx + ratio
                k = k+1
            end do
            ry = ry + ratio
        end do

        v(:,1) = [(rn(),i=1,N)]
        v(:,2) = [(rn(),i=1,N)]

        call accel(r,a)
        
    end subroutine init_md

    subroutine verlet
        implicit none        
        r = r + v*dt + 0.5d0 * a*dt2
        call pbc(r)
        v = v + 0.5d0 * a*dt

        call accel(r,a)
        v = v + 0.5d0 * a*dt
    end subroutine verlet
end module md

module measure_mod
    use md
    implicit none
    real(8)::E, Ek, Ep, vc(2), Te, P, nleft, virial
    !   E: energy
    !   Ek: kinetic energy
    !   Ep: potential energy
    !   vc: center moment, sum(v)/N
    !   T: temperature
    !   P: pressure

contains
    subroutine measure
        implicit none
        integer::i,j
        real(8)::area,dr(2)
        
        vc = sum(v)/N

        Ek = 0.5d0*sum(v*v)
        
        Ep = 0.d0
        do i = 1, N-1
            do j = i+1, N
                dr = disp(r(i,:),r(j,:))
                Ep = Ep + LJ_potential(dr)
                virial = virial + sum(dr * force(dr))
            end do
        end do

        E = Ep + Ek

        Te = Ek/(N-1)

        area = product(L)
        P = N*Te/area + virial/2/area

        nleft =1 - sum(int(r(:,1)/L(1)*2))/dble(N)

    end subroutine measure
    
end module measure_mod

program main
    use measure_mod
    implicit none
    real(8)::t

    call init_md
    open(10,file=trim(fn),status='replace')
    t = 0.d0
    do while(t<=tend)
        call verlet
        call measure
        write(10,'(9f20.8)')t, E, Ek, Ep, Te, P, vc, nleft
        t = t + dt
    end do

end program main