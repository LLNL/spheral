      subroutine uppers(hoc,next,x,up,delta,length,periodic,
     >  maxxy,maxxz)
c     
      implicit none
c     
      integer maxxy,maxxz
      integer hoc,next(maxxy),x(maxxz),up(maxxz),delta,length
      integer point,point_d,last
      logical periodic
c     
      point=hoc
      do while(point .gt. 0)
        up(point)=-1
        point=next(point)
      end do
c     
      if(periodic) then
        point=hoc
        last=-1
        if(point .gt. 0) up(point)=-1
        do while(point .gt. 0)
          point_d=next(point)
          if(point_d .gt. 0) then
            if(mod(x(point_d)+delta+length,length) .eq. 
     >        mod(x(point)+length,length)) up(point_d)=point
          end if
          last=point
c     write(59,*)' d ',point
          point=next(point)
        end do
        if(hoc .gt. 0) then
          if(mod(x(hoc)+delta+length,length) .eq. 
     >      mod(x(last)+length,length)) up(hoc)=last
        end if
      else
        point=hoc
        do while(point .gt. 0)
          point_d=next(point)
          if(point_d .gt. 0) then
            if(x(point_d)+delta .eq. x(point)) up(point_d)=point
          end if
c     write(59,*)' e ',point
          point=next(point)
        end do
      end if
      return
      end
      
