      subroutine backwards(group,hoc_a,next_a,hoc_b,next_b)
c     
      implicit none
c     
      integer hoc_a(*),hoc_b,next_a(*),next_b(*),point,group
c     
      hoc_b=-1
      point=hoc_a(group)
      do while(point .gt. 0)
       next_b(point)=hoc_b
       hoc_b=point
       point=next_a(point)
      end do
      return
      end
