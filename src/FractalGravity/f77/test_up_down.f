      subroutine test_up_down(hoc,group,level,next,pointer,
     >     up_x,up_y,up_z,
     >     down_x,down_y,down_z,x,y,z,points_maxx)
c
      implicit none
c
      integer points_maxx
      integer hoc,next(points_maxx),pointer(points_maxx)
      integer up_x(points_maxx),up_y(points_maxx)
      integer up_z(points_maxx),down_x(points_maxx),down_y(points_maxx)
      integer down_z(points_maxx)
      integer x(points_maxx),y(points_maxx),z(points_maxx)
      integer point,pp,pp_x,pp_y,pp_z,group,level
c
      logical p_x_up,p_y_up,p_z_up,p_x_down,p_y_down,p_z_down
      logical must_x_up,must_y_up,must_z_up
      logical must_x_down,must_y_down,must_z_down,bad
      logical bad_up_x,bad_up_y,bad_up_z
      logical bad_down_x,bad_down_y,bad_down_z
c
      point=hoc
      do while(point .gt. 0) 
         p_x_up=up_x(point) .lt. 0
         p_y_up=up_y(point) .lt. 0
         p_z_up=up_z(point) .lt. 0
         p_x_down=down_x(point) .lt. 0
         p_y_down=down_y(point) .lt. 0
         p_z_down= down_z(point) .lt. 0
         pp=-pointer(point)
         pp_x=mod(pp-1,3)+1
         pp_y=mod((pp-1)/3,3)+1
         pp_z=(pp-1)/9+1
         must_x_up=pp_x .le. 2
         must_y_up=pp_y .le. 2
         must_z_up=pp_z .le. 2
         must_x_down=pp_x .ge. 2
         must_y_down=pp_y .ge. 2
         must_z_down=pp_z .ge. 2
         bad_up_x=must_x_up .and. p_x_up
         bad_up_y=must_y_up .and. p_y_up
         bad_up_z=must_z_up .and. p_z_up
         bad_down_x=must_x_down .and. p_x_down
         bad_down_y=must_y_down .and. p_y_down
         bad_down_z=must_z_down .and. p_z_down
         bad=bad_up_x .or. bad_up_y .or. bad_up_z .or. 
     >     bad_down_x .or. bad_down_y .or. bad_down_z
         if(bad)write(94,94)group,level,point,x(point),y(point),
     >        z(point),pointer(point),up_x(point),up_y(point),
     >     up_z(point),down_x(point),down_y(point),down_z(point),
     >     bad_up_x,bad_up_y,bad_up_z,bad_down_x,bad_down_y,bad_down_z
 94      format(13i8,1x,6l1)
         point=next(point)
      end do
c
      return
      end
