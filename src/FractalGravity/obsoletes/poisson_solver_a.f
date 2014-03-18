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
      integer tweaks,n_xy,n_xz,n_yz,delta,grid_multiply
      integer pos_point_x(points_maxx),pos_point_y(points_maxx)
      integer pos_point_z(points_maxx)
c     
      logical inside(points_maxx),debug,periodic
c     
      real potential_point(points_maxx),density(points_maxx)
      real grav_const,eps,rjac,pi,t_1,t_2
c     
      data pi/3.141592653589793/
      if(debug)print*,'here in poisson solver',points_maxx
c     
      call cpu_time(t_1)
c
      grid_multiply=grid_length*2**level_max
      grav_const=4.0*pi/(float(grid_length*2**level(group)))**2
c     
      eps=1.0e-5
      maxits=1000
c     
      hoc_hoc_up_x=-1
c     
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
c     
c     point=hoc_hoc_up_x
c     do while(point .gt. 0)
c     write(93,*)'hoc_hoc ',group,point,inside(point)
c     point=next_hoc_up_x(point)
c     end do
c     
c      width=(float(n_tot))**(1.0/3.0)
c      rjac=cos(3.1415926535/width)     
c
      delta=2**(level_max-level(group))
c      call better_order(pos_point_z,pos_point_y,points_maxx,
c     >  hoc_hoc_up_x,next_hoc_up_x,grid_multiply,delta,periodic)
c
      rjac=(
     >  cos(3.1415926535*float(n_yz)/float(n_tot))+
     >  cos(3.1415926535*float(n_xz)/float(n_tot))+
     >  cos(3.1415926535*float(n_xy)/float(n_tot)))/3.0
c
      call sor(group,hoc_points,next_points,
     >  hoc_hoc_up_x,next_hoc_up_x,
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
