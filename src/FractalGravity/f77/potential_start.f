      subroutine potential_start(group,hoc_points,next_points,
     >  potential_point,point_up_x,point_up_y,point_up_z,
     >  point_down_x,point_down_y,point_down_z,
     >  point_pointer,grid_length,level,inside,
     >  debug,groups_maxx,points_maxx,particles_maxx,tweaks)
c     
      implicit none
c     
      integer groups_maxx,points_maxx,particles_maxx
      logical inside(points_maxx),debug
c     
      real potential_point(points_maxx)
c     
      integer point,hoc_points(groups_maxx),group,grid_length
      integer level(groups_maxx)
      integer point_up_x(points_maxx),point_up_y(points_maxx)
      integer point_up_z(points_maxx)
      integer p_up_x,p_down_x,tweaks
      integer point_down_x(points_maxx),point_down_y(points_maxx)
      integer point_down_z(points_maxx)
      integer p_up_y,p_down_y,p_up_z,p_down_z
      integer point_pointer(points_maxx),next_points(points_maxx),p_p
c     
      if(debug)write(41,*)'enter potential start ',group,level(group)

      point=hoc_points(group)
      do while(point .gt. 0)
       if(point_pointer(point) .gt. 0) then
        potential_point(point)=potential_point(point_pointer(point))
c     
        if(debug) write(41,*)'a ',group,point,point_pointer(point),
     >    potential_point(point)
       end if
       point=next_points(point)
      end do
c     
      point=hoc_points(group)
      do while(point .gt. 0)
       p_up_x=point_up_x(point)
       p_up_y=point_up_y(point)
       p_up_z=point_up_z(point)
       p_down_x=point_down_x(point)
       p_down_y=point_down_y(point)
       p_down_z=point_down_z(point)
       p_p=-point_pointer(point)
c     
       if(p_p .eq. 2 .or. p_p .eq. 8 .or.
     >   p_p .eq. 20 .or. p_p .eq. 26) then
        potential_point(point)=0.5*(
     >    potential_point(p_up_x)+potential_point(p_down_x))
c     
       else if(p_p .eq. 4 .or. p_p .eq. 6 .or. 
     >    p_p .eq. 22 .or. p_p .eq. 24) then
        potential_point(point)=0.5*(
     >    potential_point(p_up_y)+potential_point(p_down_y))
c     
       else if(p_p .eq. 10 .or. p_p .eq. 12 .or. 
     >    p_p .eq. 16 .or. p_p .eq. 18) then
        potential_point(point)=0.5*(
     >    potential_point(p_up_z)+potential_point(p_down_z))
c     
       else if(p_p .eq. 5 .or. p_p .eq. 23) then
        potential_point(point)=0.25*(
     >    potential_point(point_down_x(p_down_y))+
     >    potential_point(point_up_x(p_down_y))+
     >    potential_point(point_down_x(p_up_y))+
     >    potential_point(point_up_x(p_up_y)))
       else if(p_p .eq. 11 .or. p_p .eq. 17) then
        potential_point(point)=0.25*(
     >    potential_point(point_down_x(p_down_z))+
     >    potential_point(point_up_x(p_down_z))+
     >    potential_point(point_down_x(p_up_z))+
     >    potential_point(point_up_x(p_up_z)))
       else if(p_p .eq. 13 .or. p_p .eq. 15) then
        potential_point(point)=0.25*(
     >    potential_point(point_down_y(p_down_z))+
     >    potential_point(point_up_y(p_down_z))+
     >    potential_point(point_down_y(p_up_z))+
     >    potential_point(point_up_y(p_up_z)))
       else if(p_p .eq. 14) then
        potential_point(point)=0.125*(
     >    potential_point(
     >    point_down_z(point_down_x(p_down_y)))+
     >    potential_point(
     >    point_down_z(point_up_x(p_down_y)))+
     >    potential_point(
     >    point_down_z(point_down_x(p_up_y)))+
     >    potential_point(
     >    point_down_z(point_up_x(p_up_y)))+
     >    potential_point(
     >    point_up_z(point_down_x(p_down_y)))+
     >    potential_point(
     >    point_up_z(point_up_x(p_down_y)))+
     >    potential_point(
     >    point_up_z(point_down_x(p_up_y)))+
     >    potential_point(
     >    point_up_z(point_up_x(p_up_y))))
       else if(p_p .ge. 0) then
        stop 'potential not with you at point'
       end if
c     
       point=next_points(point)
      end do
      if(debug) then
       point=hoc_points(group)
       do while(point .gt. 0)
        write(41,*)'pot start ',group,point,point_pointer(point),
     >    potential_point(point)
        point=next_points(point)
       end do
      end if
      return
      end
