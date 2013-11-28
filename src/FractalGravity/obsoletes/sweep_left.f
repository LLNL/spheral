      subroutine sweep_left(hoc_hoc_up,next_hoc_up,n_tot,
     >  group,hoc_points,next_points,point_up,inside,
     >  groups_maxx,points_maxx)
c     
      implicit none
c     
      include 'maxx.inc'
c     
      integer groups_maxx,points_maxx
      integer hoc_hoc_up,next_hoc_up(maxx),group,n_tot
      integer hoc_points(groups_maxx),next_points(points_maxx)
      integer point_up(points_maxx),point,p_up
c     
      logical inside(points_maxx)
c     
      hoc_hoc_up=-1
      n_tot=0
      point=hoc_points(group)
      do while(point .gt. 0)
        n_tot=n_tot+1
        if(.not. inside(point)) then
          p_up=point_up(point)
          if(p_up .gt. 0) then
            if(inside(p_up)) then
c              write(82,*)p_up
              next_hoc_up(p_up)=hoc_hoc_up
              hoc_hoc_up=p_up
            end if
          end if
        end if
        point=next_points(point)
      end do
c     
      return
      end
