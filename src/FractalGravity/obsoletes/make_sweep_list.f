      subroutine make_sweep_list(hoc_hoc_up,next_hoc_up,
     >  point_up,hoc,update,inside,points_maxx)
c     
      implicit none
c
      include 'maxx.inc'
c
      integer points_maxx
      integer hoc_hoc_up,next_hoc_up(maxx),point_up(points_maxx)
      integer point,hoc(2),update(maxx),ipass,point_left
c
      logical inside(points_maxx),it_is_inside
c
      do ipass=1,2
        hoc(ipass)=-1
        point_left=hoc_hoc_up
        do while(point_left .gt. 0)
          point=point_left
          if(ipass .eq. 2) point=point_up(point)
          do while(it_is_inside(point,inside))
            update(point)=hoc(ipass)
            hoc(ipass)=point
            point=point_up(point)
            if(point .gt. 0) point=point_up(point)
          end do
          point_left=next_hoc_up(point_left)
        end do
      end do
c     
      return
      end
