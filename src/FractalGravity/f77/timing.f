      subroutine timing(what,which)
c     
      implicit none
c     
      integer what,which,n
      integer steps
      real time(30),time_0(30),time_1(30),time_2(30),t_30
      character*22 text(30)
      save time_1,time_2
      data time/30*0.0/
      data time_0/30*0.0/
      data steps/0/
      data (text(n),n=1,20)/
     >  ' offset              ',
     >  ' negatives           ',
     >  ' tree start          ',
     >  ' density 0           ',
     >  ' periodic solver     ',
     >  ' isolated solver     ',
     >  ' force at point 0    ',
     >  ' high pairs          ',
     >  ' equivalence class   ',
     >  ' high groups         ',
     >  ' daughter group      ',
     >  ' highest level group ',
     >  ' force at particle 0  ',
     >  ' find chains         ',
     >  ' densities           ',
     >  ' potential start     ',
     >  ' poisson solver      ',
     >  ' force at point tree ',
     >  ' force at parts tree  ',
     >  ' dump everyhing    '
     >  /
c     
      data text(30)/
     >  ' total time          '/
c     
      if(what .eq. -1) then
c     call second(time_1(which))
         call cpu_time(time_1(which))
         return
      else if(what .eq. 1) then
c     call second(time_2(which))
         call cpu_time(time_2(which))
         time(which)=time(which)+time_2(which)-time_1(which)
         return
      else if(what .eq. 0) then
         steps=steps+1
         t_30=(time(30)-time_0(30))+1.0e-30
         write(54,15)
         write(54,114)steps,text(30),time(30),t_30
         write(54,15)
         do n=1,20
            write(54,14)text(n),time(n),
     >           100.0*time(n)/(time(30)+1.0e-30),
     >           time(n)-time_0(n),100.0*(time(n)-time_0(n))/t_30
 114        format(i5,a22,6f10.4)
 14         format(a22,6f10.4)
 15         format(a1)
         end do
c     
         do n=1,30
            time_0(n)=time(n)
         end do
      else
         stop 'die timing'
      end if
      return
      end
