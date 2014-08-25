      logical function it_is_inside(point,inside)
c
      implicit none
c
      logical inside(*)
      integer point
c
      it_is_inside=point .gt. 0
      if(it_is_inside) it_is_inside=inside(point)
c
      return
      end
