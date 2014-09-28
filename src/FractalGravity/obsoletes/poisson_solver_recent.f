      subroutine poisson_solver(group,hoc_points,next_points,
     >  pos_point_x,pos_point_y,pos_point_z,
     >  point_up_x,point_up_y,point_up_z,
     >  point_down_x,point_down_y,point_down_z,periodic,
     >  inside,potential_point,density,level,level_max,grid_length,
     >  debug,groups_maxx,points_maxx,particles_maxx,tweaks)
c     
      implicit none
c     
      integer groups_maxx,points_maxx,particles_maxx,maxits,tweaks
c     
      include 'maxx.inc'
c     
      integer group,hoc_points(groups_maxx),next_points(points_maxx)
      integer level(groups_maxx),level_max,grid_length,hoc_a(2),hoc_b(2)
      integer point_up_x(points_maxx),point_up_y(points_maxx)
      integer point_up_z(points_maxx)
      integer point_down_x(points_maxx),point_down_y(points_maxx)
      integer point_down_z(points_maxx),update_a(maxx),update_b(maxx)
      integer hoc_hoc_up,next_hoc_up(maxx),n_tot,coord
      integer hoc_b_tmp,grid_multiply,delta
      integer pos_point_x(points_maxx),pos_point_y(points_maxx)
      integer pos_point_z(points_maxx)
c     
      logical periodic,inside(points_maxx),wrapped,debug,testing
c     
      real potential_point(points_maxx),density(points_maxx)
      real grav_const,eps,rjac,width,pi,t_1,t_2
c     
      data pi/3.141592653589793/
c     print*,'here in poisson solver',points_maxx
c     
      call cpu_time(t_1)
c
      grav_const=4.0*pi/(float(grid_length*2**level(group)))**2
c     
      grid_multiply=grid_length*2**level_max
      delta=2**(level_max-level(group))
c     
      eps=1.0e-5
      maxits=500
      coord=1
      wrapped=.false.
      testing=.false.
c     
c     write(82,*)'group ',group
      call sweep_left(hoc_hoc_up,next_hoc_up,n_tot,
     >  group,hoc_points,next_points,point_up_x,inside,
     >  groups_maxx,points_maxx)
c     
      coord=1
c     
      if(periodic) then
c     
        call test_wrap(group,hoc_points,next_points,inside,
     >    hoc_hoc_up,next_hoc_up,point_up_x,wrapped,
     >    groups_maxx,points_maxx)
c     
        if(wrapped)then
c     
          call sweep_left(hoc_hoc_up,next_hoc_up,n_tot,
     >      group,hoc_points,next_points,point_up_y,inside,
     >      groups_maxx,points_maxx)
c     
          call test_wrap(group,hoc_points,next_points,inside,
     >      hoc_hoc_up,next_hoc_up,point_up_y,wrapped,
     >      groups_maxx,points_maxx)
c     
          coord=2
        end if
c     
        if(wrapped)then
          call sweep_left(hoc_hoc_up,next_hoc_up,n_tot,
     >      group,hoc_points,next_points,point_up_z,inside,
     >      groups_maxx,points_maxx)
c     
          call test_wrap(group,hoc_points,next_points,inside,
     >      hoc_hoc_up,next_hoc_up,point_up_z,wrapped,
     >      groups_maxx,points_maxx)
c     
          coord=3
        end if
        if(wrapped) then
          print*,'wrapped ',group
          stop 'die wrapped from poisson'
        end if
        if(coord .ne. 1) print*,'wrapped ',group,coord
      end if
c     
      if(testing) then
        if(n_tot .ge. 1000) then
          if(coord .eq. 1) then
            call better_order(pos_point_y,pos_point_z,points_maxx,
     >        hoc_hoc_up,next_hoc_up,grid_multiply,delta,periodic)
          else if(coord .eq. 2) then
            call better_order(pos_point_x,pos_point_z,points_maxx,
     >        hoc_hoc_up,next_hoc_up,grid_multiply,delta,periodic)
          else
            call better_order(pos_point_x,pos_point_y,points_maxx,
     >        hoc_hoc_up,next_hoc_up,grid_multiply,delta,periodic)
          end if
        end if
c     
      end if
      if(coord .eq. 1) then
        call make_sweep_list(hoc_hoc_up,next_hoc_up,
     >    point_up_x,hoc_a,update_a,inside,points_maxx)
      else if(coord .eq. 2) then
        call make_sweep_list(hoc_hoc_up,next_hoc_up,
     >    point_up_y,hoc_a,update_a,inside,points_maxx)
      else
        call make_sweep_list(hoc_hoc_up,next_hoc_up,
     >    point_up_z,hoc_a,update_a,inside,points_maxx)
      end if
c     
c     
      call backwards(1,hoc_a,update_a,hoc_b_tmp,update_b)
      hoc_b(1)=hoc_b_tmp
      call backwards(2,hoc_a,update_a,hoc_b_tmp,update_b)
      hoc_b(2)=hoc_b_tmp
c     
      width=(float(n_tot))**(1.0/3.0)
      rjac=cos(3.1415926535/width)
c     
      call sor(group,hoc_a,update_a,hoc_b,update_b,
     >  potential_point,density,point_up_x,point_up_y,point_up_z,
     >  point_down_x,point_down_y,point_down_z,grav_const,
     >  rjac,eps,debug,maxits,points_maxx)
c     
      call cpu_time(t_2)
      write(93,*)'n_tot ',group,level(group),n_tot,coord,t_2-t_1
c     
      return
      end
