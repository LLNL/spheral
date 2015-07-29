      subroutine negatives_points(points_maxx,high_points_maxx,
     >  points_maxx_group,next_points,hoc_particles,pos_point_x,
     >  pos_point_y,pos_point_z,point_up_x,point_down_x,point_up_y,
     >  point_down_y,point_up_z,point_down_z,inside,
     >  particles_at_point,density,
     >  potential_point,force_point_x,force_point_y,force_point_z,
     >  seq_to_list,
     >  list_high_1,list_high_2,group_tmp,next_high_points,tweaks)
c     
      implicit none
c
      include 'maxx.inc'
c     
      integer points_maxx,high_points_maxx,point,tweaks
      integer points_maxx_group,next_points(points_maxx)
c
      real density(points_maxx)
      real potential_point(points_maxx),force_point_x(points_maxx)
      real force_point_y(points_maxx),force_point_z(points_maxx)
c     
      integer hoc_particles(points_maxx)
      integer pos_point_x(points_maxx),pos_point_y(points_maxx)
      integer pos_point_z(points_maxx)
      integer point_up_x(points_maxx)
      integer point_down_x(points_maxx),point_up_y(points_maxx)
      integer point_down_y(points_maxx)
      integer point_up_z(points_maxx)
      integer point_down_z(points_maxx)
      integer seq_to_list(high_points_max_group)
      integer list_high_1(3*high_points_max_group)
      integer list_high_2(3*high_points_max_group)
      integer group_tmp(high_points_max_group)
      integer next_high_points(high_points_maxx)
      integer particles_at_point(points_maxx)
c     
      logical inside(points_maxx)
c
      do point=1,points_maxx
       next_points(point)=-point
       hoc_particles(point)=-point
       pos_point_x(point)=-point
       pos_point_y(point)=-point
       pos_point_z(point)=-point
       point_up_x(point)=-point
       point_up_y(point)=-point
       point_up_z(point)=-point
       point_down_x(point)=-point
       point_down_y(point)=-point
       point_down_z(point)=-point
       particles_at_point(point)=-point
       force_point_x(point)=-float(point)*1.0e30
       force_point_y(point)=-float(point)*1.0e30
       force_point_z(point)=-float(point)*1.0e30
       density(point)=0.0
       potential_point(point)=-float(point)*1.0e30
       inside(point)=.false.
      end do
c     
      do point=1,high_points_maxx
       next_high_points(point)=-point
      end do
c     
      do point=1,high_points_max_group
         seq_to_list(point)=-point
         group_tmp(point)=-point
         list_high_1(point)=-point
         list_high_2(point)=-point
      end do
c
      do point=high_points_max_group+1,3*high_points_max_group
         list_high_1(point)=-point
         list_high_2(point)=-point
      end do
c
      return
      end
