      subroutine force_at_point(group,hoc_points,next_points,
     >  force_point_x,force_point_y,force_point_z,
     >  potential_point,point_up_x,point_up_y,point_up_z,
     >  point_down_x,point_down_y,point_down_z,
     >  point_pointer,grid_length,level,inside,
     >  debug,groups_maxx,points_maxx,particles_maxx,tweaks)
c     
      implicit none
c     
      integer groups_maxx,points_maxx,particles_maxx,tweaks
      logical inside(points_maxx),debug
c     
      real conv,force_point_x(points_maxx),force_point_y(points_maxx)
      real force_point_z(points_maxx)
      real potential_point(points_maxx)
c     
      integer point,hoc_points(groups_maxx),group,grid_length
      integer level(groups_maxx)
      integer point_up_x(points_maxx),point_up_y(points_maxx)
      integer point_up_z(points_maxx)
      integer p_up_x,p_down_x,p_up_y,p_down_y,p_up_z,p_down_z,p_p
      integer point_down_x(points_maxx),point_down_y(points_maxx)
      integer point_down_z(points_maxx)
      integer point_pointer(points_maxx),next_points(points_maxx)
      integer p_p_up,p_p_down,p_p_up_up,p_p_down_down
      integer p_p_down_up,p_p_up_down
c     
      conv=float(grid_length)*2.0**(level(group)-1)
      if(debug)
     >  write(34,*)'enter force at point ',group,level(group),conv
c     
      point=hoc_points(group)
      do while(point .gt. 0)
       p_p=-point_pointer(point)
       p_up_x=point_up_x(point)
       p_up_y=point_up_y(point)
       p_up_z=point_up_z(point)
       p_down_x=point_down_x(point)
       p_down_y=point_down_y(point)
       p_down_z=point_down_z(point)
c     
c     if(debug) write(34,*)'force ',point,p_up_x,p_up_y,p_up_z,
c     >        p_down_x,p_down_y,p_down_z,p_p,inside(point)
c     
       if(inside(point)) then
        force_point_x(point)=(potential_point(p_down_x)-
     >    potential_point(p_up_x))*conv
        force_point_y(point)=(potential_point(p_down_y)-
     >    potential_point(p_up_y))*conv
        force_point_z(point)=(potential_point(p_down_z)-
     >    potential_point(p_up_z))*conv
c     
       else if(p_p .lt. 0) then
        force_point_x(point)=force_point_x(-p_p)
        force_point_y(point)=force_point_y(-p_p)
        force_point_z(point)=force_point_z(-p_p)
c     
       else if(p_p .eq. 2 .or. p_p .eq. 8 .or.
     >    p_p .eq. 20 .or. p_p .eq. 26) then
        p_p_up=point_pointer(p_up_x)
        p_p_down=point_pointer(p_down_x)
c     write(34,*)'what ',p_up_x,p_p_up,p_down_x,p_p_down
        force_point_x(point)=0.5*(
     >    force_point_x(p_p_up)+force_point_x(p_p_down))
        force_point_y(point)=0.5*(
     >    force_point_y(p_p_up)+force_point_y(p_p_down))
        force_point_z(point)=0.5*(
     >    force_point_z(p_p_up)+force_point_z(p_p_down))
c     
       else if(p_p .eq. 4 .or. p_p .eq. 6 .or. 
     >    p_p .eq. 22 .or. p_p .eq. 24) then
        p_p_up=point_pointer(p_up_y)
        p_p_down=point_pointer(p_down_y)
        force_point_x(point)=0.5*(
     >    force_point_x(p_p_up)+force_point_x(p_p_down))
        force_point_y(point)=0.5*(
     >    force_point_y(p_p_up)+force_point_y(p_p_down))
        force_point_z(point)=0.5*(
     >    force_point_z(p_p_up)+force_point_z(p_p_down))
c     
       else if(p_p .eq. 10 .or. p_p .eq. 12 .or. 
     >    p_p .eq. 16 .or. p_p .eq. 18) then
        p_p_up=point_pointer(p_up_z)
        p_p_down=point_pointer(p_down_z)
        force_point_x(point)=0.5*(
     >    force_point_x(p_p_up)+force_point_x(p_p_down))
        force_point_y(point)=0.5*(
     >    force_point_y(p_p_up)+force_point_y(p_p_down))
        force_point_z(point)=0.5*(
     >    force_point_z(p_p_up)+force_point_z(p_p_down))
c     
       else if(p_p .eq. 5 .or. p_p .eq. 23) then
        p_p_up_up=point_pointer(point_up_x(p_up_y))
        p_p_down_up=point_pointer(point_down_x(p_up_y))
        p_p_down_down=point_pointer(point_down_x(p_down_y))
        p_p_up_down=point_pointer(point_up_x(p_down_y))
        force_point_x(point)=0.25*(
     >    force_point_x(p_p_down_down)+force_point_x(p_p_up_down)+
     >    force_point_x(p_p_down_up)+force_point_x(p_p_up_up))
        force_point_y(point)=0.25*(
     >    force_point_y(p_p_down_down)+force_point_y(p_p_up_down)+
     >    force_point_y(p_p_down_up)+force_point_y(p_p_up_up))
        force_point_z(point)=0.25*(
     >    force_point_z(p_p_down_down)+force_point_z(p_p_up_down)+
     >    force_point_z(p_p_down_up)+force_point_z(p_p_up_up))
       else if(p_p .eq. 11 .or. p_p .eq. 17) then
        p_p_up_up=point_pointer(point_up_x(p_up_z))
        p_p_down_up=point_pointer(point_down_x(p_up_z))
        p_p_down_down=point_pointer(point_down_x(p_down_z))
        p_p_up_down=point_pointer(point_up_x(p_down_z))
        force_point_x(point)=0.25*(
     >    force_point_x(p_p_down_down)+force_point_x(p_p_up_down)+
     >    force_point_x(p_p_down_up)+force_point_x(p_p_up_up))
        force_point_y(point)=0.25*(
     >    force_point_y(p_p_down_down)+force_point_y(p_p_up_down)+
     >    force_point_y(p_p_down_up)+force_point_y(p_p_up_up))
        force_point_z(point)=0.25*(
     >    force_point_z(p_p_down_down)+force_point_z(p_p_up_down)+
     >    force_point_z(p_p_down_up)+force_point_z(p_p_up_up))
       else if(p_p .eq. 13 .or. p_p .eq. 15) then
        p_p_up_up=point_pointer(point_up_y(p_up_z))
        p_p_down_up=point_pointer(point_down_y(p_up_z))
        p_p_down_down=point_pointer(point_down_y(p_down_z))
        p_p_up_down=point_pointer(point_up_y(p_down_z))
        force_point_x(point)=0.25*(
     >    force_point_x(p_p_down_down)+force_point_x(p_p_up_down)+
     >    force_point_x(p_p_down_up)+force_point_x(p_p_up_up))
        force_point_y(point)=0.25*(
     >    force_point_y(p_p_down_down)+force_point_y(p_p_up_down)+
     >    force_point_y(p_p_down_up)+force_point_y(p_p_up_up))
        force_point_z(point)=0.25*(
     >    force_point_z(p_p_down_down)+force_point_z(p_p_up_down)+
     >    force_point_z(p_p_down_up)+force_point_z(p_p_up_up))
       else 
        stop 'force not with you at this point'
       end if
       point=next_points(point)
      end do

      return
      end
