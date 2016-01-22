      subroutine better_order_1(pos_b,points_maxx,
     >  hoc_hoc_up,next_hoc_up,grid_multiply,delta,periodic)
c     
      implicit none
c
      include 'maxx.inc'
c
      integer points_maxx
      integer pos_b(points_maxx)
      integer hoc_hoc_up,next_hoc_up(points_maxx)
      integer grid_multiply,delta,b_min,b_max,hoc_left,total_b
      integer n,hoc_b(maxx),next_b(maxx),point
c
      logical periodic
c     
      if(periodic) return
c     
      b_min=grid_multiply
      b_max=-b_min
      hoc_left=hoc_hoc_up
      do while(hoc_left .gt. 0)
        b_min=min(b_min,pos_b(hoc_left))
        b_max=max(b_max,pos_b(hoc_left))
        hoc_left=next_hoc_up(hoc_left)
      end do
c
      total_b=(b_max-b_min)/delta+1
      do n=1,total_b
        hoc_b(n)=-1
      end do
c
      hoc_left=hoc_hoc_up
      do while(hoc_left .gt. 0)
        n=(pos_b(hoc_left)-b_min)/delta+1
        next_b(hoc_left)=hoc_b(n)
        hoc_b(n)=hoc_left
        hoc_left=next_hoc_up(hoc_left)
      end do
c
      hoc_hoc_up=-1
      do n=1,total_b
        point=hoc_b(n)
        do while(point .gt. 0)
          next_hoc_up(point)=hoc_hoc_up
          hoc_hoc_up=point
          point=next_b(point)
        end do
      end do
c     
      return
      end
