! metropolis and wolff cluster MC of 2D Ising model
! version 2016, by Wenan Guo
      module system
      implicit none
      integer, parameter:: maxl=1024,maxsize=maxl*maxl
      real(8):: ppr(-4:4),ppm(-4:4)
      real(8):: bp
      integer :: isp(maxsize),nbor(4,maxsize),n1,nsq,msign
      endmodule system
      module data
      real(8):: arm, asm, ave, ase
      endmodule data

      program ising
      use system; use data
      implicit none
      integer:: iseed,msample,nbin,neq,nh,nm,nw,i,ib,ip
      real(8):: temp
      save
 10   continue                                                         
      read *,n1,neq,iseed,msample,nbin,nh,nm,nw,temp
      if ((n1.eq.0).or.(n1.gt.maxl).or.(iseed.eq.0)) stop
 12   continue     
      print *, 'simulate temperature=',temp, ' size=',n1
! initialize, equilibrate
      call initial(temp,iseed)
      do i=1,neq
        call simul(nh,nm,nw)
      enddo
      do ib=1,nbin
        call cleardata
!       simulate and sample
        do ip=1,msample                                                   
          call simul(nh,nm,nw)                          
          call sample
        enddo                                                          
        call writeres(msample,temp)
      enddo
      read *,nh,nm,nw,temp
      if (temp.eq.0.0d0) goto 10       ! next run
      goto 12                          ! next subrun
      end
      subroutine writeres(msample,temp)
      use data; use system; implicit none
      integer:: msample
      real(8):: c,q,temp
! normalize, analyze, and print output
      arm=arm/msample
      asm=asm/msample
      ave=ave/msample
      ase=ase/msample
      c=(ase-ave**2)/temp**2*n1**2
      q=asm/arm**2
      open(7, file='bin.dat',access='append')
      write(7, '(i5,6f15.8)')n1, temp, arm, asm, ave, c, q
      close(7)
      endsubroutine writeres

      subroutine cleardata
      use data
      arm=0              ! average relative magnetization
      asm=0              ! average square magnetization
      ave=0              ! average relative energy
      ase=0              ! average square energy
      endsubroutine cleardata

      subroutine initial(temp,iseed)
!     initialize for Ising mc                                           
      use system
      implicit real*8(a-h,o-z)                                          
      save
      nsq=n1*n1                                                        
      print * , 'temp=',temp
! define neighbors of each spin
      do 210 ispin=1,nsq                                                
      iy=(ispin-1)/n1+1                                                
      ix=ispin-(iy-1)*n1                                     
      ixp=ix+1-(ix/n1)*n1                                               
      iyp=iy+1-(iy/n1)*n1                                               
      ixm=ix-1+((n1-ix+1)/n1)*n1                                        
      iym=iy-1+((n1-iy+1)/n1)*n1                                        
      nbor(1,ispin)=(iy -1)*n1+ixm 
      nbor(4,ispin)=(iy -1)*n1+ixp 
      nbor(2,ispin)=(iym-1)*n1+ix  
      nbor(3,ispin)=(iyp-1)*n1+ix
210   continue                                                          
! fill lookup table for local updates and define bond probabilities
      bfm=dexp(1.0d0/temp) 
! bond probability for wolff cluster
      bp=1.0d0-1.0d0/(bfm*bfm)
! heat bath prob.
      do 220 i=0,4                                                   
      ppr(i)=1.d0/(1.d0+dexp(-i*2.d0/temp))
      ppr(-i)=1.d0/(1.d0+dexp(i*2.d0/temp))
220   continue                                                          
      do i=0,4
      ppm(i)=dexp(-i*2.d0/temp)
      ppm(-i)=1.d0
      enddo
! initialize random generator
      call ransi(iabs(iseed)+1)
! initialize spin config 
      is=1
      do 112 ns=1,nsq
      isp(ns)=is
112   is=-is
      return                                                            
      end                                                               

      subroutine simul(nhsteps,nmsteps,nwsteps)                
      implicit none
      integer::nmsteps,nwsteps,nhsteps
!     execute mc steps
      call mch(nhsteps)
      call mcm(nmsteps)
      call mcw(nwsteps)
      return                                                            
      end                                                               

      subroutine mch(nsteps)                                     
! subroutine for heat bath sweeps
      use system
      implicit none
      integer::nsteps,isteps,ispin,ns,it
      real(8)::rn
      do 200 isteps=1,nsteps
      do 200 ispin=1,nsq
      ns=int(rn()*nsq)+1
      it=isp(nbor(1,ns))+isp(nbor(2,ns))+isp(nbor(3,ns))+isp(nbor(4,ns))
      if (rn().gt.ppr(it)) then 
         isp(ns)=-1
      else
         isp(ns)=1
      endif
 200  continue                                                          
      return                                                            
      end                                                               

      subroutine mcm(nsteps)                                     
! subroutine for metropolis sweeps
      use system
      implicit none
      integer::nsteps,isteps,ispin,ns,it
      real(8)::rn
      do 200 isteps=1,nsteps
      do 200 ispin=1,nsq
      ns=int(rn()*nsq)+1
      it=isp(nbor(1,ns))+isp(nbor(2,ns))+isp(nbor(3,ns))+isp(nbor(4,ns))
      if (rn().lt.ppm(it*isp(ns))) then 
         isp(ns)=-isp(ns)
      endif
 200  continue                                                          
      return                                                            
      end                                                               


      subroutine mcw(nsteps)                                     
      use system
! subroutine flips nsteps wolff clusters 
      implicit none
      integer:: nsteps,isteps,ns,icsp,nstack,js,ks
      integer:: inb
      integer:: istn(maxsize)
      real(8):: rn
      save
      do 110 isteps=1,nsteps
      ns=int(rn()*nsq)+1     
      icsp=isp(ns)
      isp(ns)=-icsp
      nstack=0
      js=ns
 104  continue  
      do 107 inb=1,4
      ks=nbor(inb,js)
      if ((isp(ks).eq.icsp).and.(rn().lt.bp)) then
         nstack=nstack+1
         istn(nstack)=ks
         isp(ks)=-icsp   
      endif
 107  continue 
      if (nstack.eq.0) goto 110
      js=istn(nstack)
      nstack=nstack-1
      goto 104
 110  continue     
      return          
      end

      subroutine sample
!     compute observables from given spin configuration                 
      use system;use data
      implicit none
      integer::m1,nbs,ns
      real(8)::enr,rem,sqm
      save
!     sample magnetization 
      m1=0
      do 20 ns=1,nsq
      m1=m1+isp(ns)
 20   continue                                                          
!     sample nearest-neighbour sum
      nbs=0
      do 30 ns=1,nsq                                                    
      nbs=nbs+isp(ns)*(isp(nbor(1,ns))+isp(nbor(2,ns)))
 30   continue                                                          
      enr=dfloat(nbs)/(n1*n1)
      rem=dabs(dfloat(m1)/(n1*n1))
!     accumulate results 
      sqm=rem*rem            ! calculate square magnetization
      arm=arm+rem            ! accumulate relative magnetization
      asm=asm+sqm            ! accumulate square magnetization
      ave=ave+enr            ! accumulate energy density
      ase=ase+enr*enr        ! accumulate square energy
      return                                                            
      end                                                               

      subroutine ransi(iseed)                                 
! initialize shift register random generator 
      implicit real*8(a-h,o-z)            
      save
      parameter (mult=32781,lenr=9689,ifdb=471)
      common/ransrb/ irs(lenr),next(lenr),ipoint,ipoinf
      k=3**18+2*iseed     
      do 100 i=1,lenr      
      k=k*mult               
      irs(i)=k+i/3+i/17       
 100  continue                 
      do 101 i=1,lenr           
 101  next(i)=i+1                
      next(lenr)=1                
      ipoint=1                     
      ipoinf=ifdb+1                
      return                         
      end

      function rn()
!     calculate random number      
      implicit real*8(a-h,o-z)            
      save
      parameter (tm32=2.d0**(-32),lenr=9689,ifdb=471)
      common/ransrb/ irs(lenr),next(lenr),ipoint,ipoinf
      irn=ieor(irs(ipoint),irs(ipoinf))   
      irs(ipoint)=irn                    
      rn=irn*tm32+0.5d0
      ipoint=next(ipoint)             
      ipoinf=next(ipoinf)            
      return                       
      end                         
