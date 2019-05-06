      program treatres
      implicit none
      integer,parameter::np=10000
      integer :: nb,i,l
      real(8) :: a(5,2),d(5),temp,t
     
      open(7,file='bin.dat')
      t=-1
      nb=10
      i=0
      a=0.d0
      do
      i=i+1 
      read(7,*,end=100) l, temp, d(1), d(2), d(3), d(4), d(5)
          if(temp.ne.t.and.t.ne.-1)then
          print *, 'temp=',temp
          print *, a(1,1), a(1,2)
          nb=i-1
          print *, 'nb=',nb
          call writeres
          a=0.d0
          i=1
      endif
          t=temp
          a(:,1)=a(:,1)+d(:)
          a(:,2)=a(:,2)+d(:)**2
      enddo
100   continue
      close(7)
      call writeres
      contains
      subroutine writeres
          a(:,1)=a(:,1)/nb
          a(:,2)=a(:,2)/nb
          a(:,2)=dsqrt(a(:,2)-a(:,1)**2)/dsqrt(nb-1.d0)
         open(8,file='m.res',access='append')
         write(8,'(i5,f10.5,2f15.8)') l, t, a(1,1),a(1,2)
         close(8)
         open(8,file='e.res',access='append')
         write(8,'(i5,f10.5,2f15.8)') l, t, a(3,1),a(3,2)
         close(8)
         open(8,file='c.res',access='append')
         write(8,'(i5,f10.5,2f15.8)') l, t, a(4,1),a(4,2)
         close(8)
         open(8,file='q.res',access='append')
         write(8,'(i5,f10.5,2f15.8)') l, t, a(5,1),a(5,2)
         close(8)
       endsubroutine writeres
      end
