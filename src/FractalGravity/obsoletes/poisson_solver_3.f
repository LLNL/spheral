      subroutine poisson_solver(group,hoc_points,next_points,
     >  pos_point_x,pos_point_y,pos_point_z,
     >  point_up_x,point_up_y,point_up_z,
     >  point_down_x,point_down_y,point_down_z,
     >  periodic,inside,potential_point,density,
     >  level,level_max,grid_length,
     >  debug,groups_maxx,points_maxx,particles_maxx,
     >  tweaks)
c     
      implicit none
c     
      integer groups_maxx,points_maxx,particles_maxx,maxits
c     
      include 'maxx.inc'
c     
      integer group,hoc_points(groups_maxx),next_points(points_maxx)
      integer level(groups_maxx),level_max,grid_length
      integer point_up_x(points_maxx),point_up_y(points_maxx)
      integer point_up_z(points_maxx)
      integer point_down_x(points_maxx),point_down_y(points_maxx)
      integer point_down_z(points_maxx)
      integer hoc_hoc_up_x,next_hoc_up_x(maxx),point,n_tot,p_up_x
      integer hoc_hoc_up_y,next_hoc_up_y(maxx),p_up_y
      integer hoc_hoc_up_z,next_hoc_up_z(maxx),p_up_z
      integer tweaks,n_xy,n_xz,n_yz
      integer pos_point_x(points_maxx),pos_point_y(points_maxx)
      integer pos_point_z(points_maxx)
c     
      logical inside(points_maxx),debug,periodic
c     
      real potential_point(points_maxx),density(points_maxx)
      real grav_const,eps,rjac,width,pi,t_1,t_2
c     
      data pi/3.141592653589793/
      if(debug)print*,'here in poisson solver',points_maxx
c     
      call cpu_time(t_1)
c
      grav_const=4.0*pi/(float(grid_length*2**level(group)))**2
c     
      eps=1.0e-5
      maxits=1000
c     
      hoc_hoc_up_x=-1
      hoc_hoc_up_y=-1
      hoc_hoc_up_z=-1
c     
      n_xy=0
      n_xz=0
      n_yz=0
      n_tot=0
      point=hoc_points(group)
      do while(point .gt. 0)
        n_tot=n_tot+1
        if(.not. inside(point)) then
          p_up_x=point_up_x(point)
          if(p_up_x .gt. 0) then
            if(inside(p_up_x)) then
              next_hoc_up_x(p_up_x)=hoc_hoc_up_x
              hoc_hoc_up_x=p_up_x
            end if
          end if
c
          p_up_y=point_up_y(point)
          if(p_up_y .gt. 0) then
            if(inside(p_up_y)) then
              next_hoc_up_y(p_up_y)=hoc_hoc_up_y
              hoc_hoc_up_y=p_up_y
            end if
          end if
c
          p_up_z=point_up_z(point)
          if(p_up_z .gt. 0) then
            if(inside(p_up_z)) then
              next_hoc_up_z(p_up_z)=hoc_hoc_up_z
              hoc_hoc_up_z=p_up_z
            end if
          end if
        end if
        if(point_down_x(point) .lt. 0) n_yz=n_yz+1
        if(point_down_y(point) .lt. 0) n_xz=n_xz+1
        if(point_down_z(point) .lt. 0) n_xy=n_xy+1
c
        point=next_points(point)
      end do
c     
      width=(float(n_tot))**(1.0/3.0)
      rjac=(
     >  cos(3.1415926535*float(n_yz)/float(n_tot))+
     >  cos(3.1415926535*float(n_xz)/float(n_tot))+
     >  cos(3.1415926535*float(n_xy)/float(n_tot)))/3.0
c
c      rjac=cos(3.1415926535/width)
c     
      call sor_3(group,hoc_points,next_points,
     >  hoc_hoc_up_x,next_hoc_up_x,
     >  hoc_hoc_up_y,next_hoc_up_y,
     >  hoc_hoc_up_z,next_hoc_up_z,
     >  potential_point,density,point_up_x,point_up_y,point_up_z,
     >  point_down_x,point_down_y,point_down_z,inside,grav_const,
     >  rjac,eps,debug,maxits,points_maxx)
c     
      call cpu_time(t_2)
      if(t_2-t_1 .gt. 0.1)
     >  write(93,*)'n_tot ',group,level(group),n_tot,t_2-t_1
c     
      return
      end
c
