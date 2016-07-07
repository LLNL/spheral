      subroutine test_wrap(group,hoc_points,next_points,inside,
     >  hoc_hoc_up,next_hoc_up,point_up,wrapped,
     >  groups_maxx,points_maxx)
c     
      implicit none
c
      integer groups_maxx,points_maxx
c
      include 'maxx.inc'
c
      integer group,hoc_points(groups_maxx),next_points(points_maxx)
      integer hoc_hoc_up,next_hoc_up(maxx),point_up(points_maxx)
      integer point,point_left
c
      logical inside(points_maxx),wrapped,touched(maxx)
c
      point=hoc_points(group)
      do while(point .gt. 0)
        if(inside(point)) touched(point)=.false.
        point=next_points(point)
      end do
c     
      point_left=hoc_hoc_up
      do while(point_left .gt. 0)
        point=point_left
        do while(inside(point))
          touched(point)=.true.
          point=point_up(point)
        end do
        point_left=next_hoc_up(point_left)
      end do
c     
      wrapped=.false.
      point=hoc_points(group)
      do while(point .gt. 0 .and. .not. wrapped)
        if(inside(point))wrapped= .not. touched(point)
        point=next_points(point)
      end do
c     
      return
      end
