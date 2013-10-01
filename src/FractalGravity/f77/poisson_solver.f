      subroutine poisson_solver(group,hoc_points,next_points,
     >  point_up_x,point_up_y,point_up_z,
     >  point_down_x,point_down_y,point_down_z,
     >  periodic,inside,potential_point,density,
     >  level,level_max,grid_length,maxits_sor,epsilon_sor,
     >  memory_value,debug,groups_maxx,points_maxx,particles_maxx,
     >  double_solver)
c     
      implicit none
c     
      integer groups_maxx,points_maxx,particles_maxx,maxits_sor
c     
      include 'maxx.inc'
c     
      integer group,hoc_points(groups_maxx),next_points(points_maxx)
      integer level(groups_maxx),level_max,grid_length
      integer point_up_x(points_maxx),point_up_y(points_maxx)
      integer point_up_z(points_maxx)
      integer point_down_x(points_maxx),point_down_y(points_maxx)
      integer point_down_z(points_maxx)
      integer point,n_tot,p_up_x,n_pointer,memory_value
      integer n_xy,n_xz,n_yz,n_inside
      integer pointer_left_x(8*grid_length_max**2)
c     
      logical inside(points_maxx),debug,periodic,double_solver
c     
      real potential_point(points_maxx),density(points_maxx)
      real grav_const,epsilon_sor,rjac,pi,t_1,t_2
c     
      if(debug)print*,'here in poisson solver',points_maxx
c     
      call cpu_time(t_1)
c
      pi=4.0*atan(1.0)
      grav_const=4.0*pi/(float(grid_length*2**level(group)))**2
c     
      n_tot=0
      n_inside=0
      n_pointer=0
      point=hoc_points(group)
      do while(point .gt. 0)
        n_tot=n_tot+1
        if(.not. inside(point)) then
          p_up_x=point_up_x(point)
          if(p_up_x .gt. 0) then
            if(inside(p_up_x)) then
               n_pointer=n_pointer+1
               pointer_left_x(n_pointer)=p_up_x
            end if
         end if
      else
         n_inside=n_inside+1
        end if
        point=next_points(point)
      end do
c     
      n_xy=0
      n_xz=0
      n_yz=0
      point=hoc_points(group)
      do while(point .gt. 0)
        if(point_down_x(point) .lt. 0) n_yz=n_yz+1
        if(point_down_y(point) .lt. 0) n_xz=n_xz+1
        if(point_down_z(point) .lt. 0) n_xy=n_xy+1
        point=next_points(point)
      end do
c     
      rjac=(
     >  cos(pi*float(n_yz)/float(n_tot))+
     >  cos(pi*float(n_xz)/float(n_tot))+
     >  cos(pi*float(n_xy)/float(n_tot)))/3.0
c
      if(double_solver) then
         call sor_d(group,hoc_points,next_points,
     >        n_pointer,pointer_left_x,
     >        potential_point,density,point_up_x,point_up_y,point_up_z,
     >        point_down_x,point_down_y,point_down_z,inside,grav_const,
     >        rjac,epsilon_sor,memory_value,debug,maxits_sor,
     >        groups_maxx,points_maxx)
      else
         call sor(group,hoc_points,next_points,
     >        n_pointer,pointer_left_x,
     >        potential_point,density,point_up_x,point_up_y,point_up_z,
     >        point_down_x,point_down_y,point_down_z,inside,grav_const,
     >        rjac,epsilon_sor,memory_value,debug,maxits_sor,
     >        groups_maxx,points_maxx)
      end if
      call cpu_time(t_2)
      if(t_2-t_1 .gt. 0.05)
     >  write(93,*)'n_tot ',group,level(group),n_inside,n_tot,t_2-t_1
c     
      return
      end
c
