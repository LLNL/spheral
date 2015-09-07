      subroutine test_it_all(grid_length,level_max,hoc_groups,
     >  next_groups,hoc_points,next_points,inside,hoc_particles,
     >  next_particles,potential_point,potential,
     >  force_point_x,force_point_y,force_point_z,particle_pointer,
     >  pos_x,pos_y,pos_z,pos_point_x,pos_point_y,pos_point_z,
     >  point_pointer,force_x,force_y,force_z,
     >  level,highest_level_group,periodic,
     >  groups_maxx,points_maxx,particles_maxx,particles_real,tweaks)
c     
      implicit none
c     
      integer groups_maxx,points_maxx,particles_maxx,particles_real
      integer grid_length,level_max,hoc_groups,group,tweaks
      integer point,hoc_points(groups_maxx),level(groups_maxx)
      integer next_points(points_maxx),next_groups(groups_maxx)
      integer particle,hoc_particles(points_maxx),part
      integer next_particles(particles_maxx)
      integer highest_level_group(particles_maxx)
      integer particle_pointer(particles_maxx)
      integer point_pointer(points_maxx)
      integer pos_point_x(points_maxx),pos_point_y(points_maxx)
      integer pos_point_z(points_maxx)
c     
      logical inside(points_maxx),periodic
c     
      real x_0,y_0,z_0,conv,d_x,d_y,d_z,d2,f_x,f_y,f_z,pot,pot_0
      real pos_x(particles_real),pos_y(particles_real)
      real pos_z(particles_real)
      real force_point_x(points_maxx),force_point_y(points_maxx)
      real force_point_z(points_maxx)
      real force_x(particles_maxx),force_y(particles_maxx)
      real force_z(particles_maxx)
      real potential_point(points_maxx)
      real potential(particles_real)
c     
      print*,'x_0,y_0,z_0 '
      read(5,*)x_0,y_0,z_0
c     
      conv=grid_length*2**level_max
c     
      group=hoc_groups
      do while(group .gt. 0)
       print*,'test it ',group,groups_maxx
       point=hoc_points(group)
       do while(point .gt. 0)
        d_x=pos_point_x(point)/conv-x_0
        d_y=pos_point_y(point)/conv-y_0
        d_z=pos_point_z(point)/conv-z_0
        d2=d_x**2+d_y**2+d_z**2+
     >    4.0**(-level(group)-1)/float(grid_length**2)
        f_x=-d_x/d2**1.5
        f_y=-d_y/d2**1.5
        f_z=-d_z/d2**1.5
        pot_0=-1.0/sqrt(d2)
        write(99,99)group,level(group),point,inside(point),
     >    float(pos_point_x(point))/
     >    float(grid_length*2**level_max),
     >    float(pos_point_y(point))/
     >    float(grid_length*2**level_max),
     >    float(pos_point_z(point))/
     >    float(grid_length*2**level_max),
     >    f_x,force_point_x(point),f_y,force_point_y(point),
     >    f_z,force_point_z(point),
     >    pot_0,potential_point(point)
 99     format(i4,i3,i8,l2,12(1pe11.3))
        if(point_pointer(point) .gt. 0) then
         if(pos_point_x(point) .ne. 
     >     pos_point_x(point_pointer(point))  .or.
     >     pos_point_y(point) 
     >     .ne. pos_point_y(point_pointer(point)))
     >     write(97,97)group,level(group),point,inside(point),
     >     point_pointer(point),
     >     pos_point_x(point),
     >     pos_point_x(point_pointer(point)),
     >     pos_point_y(point),pos_point_y(point_pointer(point))
        end if
 97     format(3i6,l2,65i6)
        point=next_points(point)
       end do
       group=next_groups(group)
      end do
c     
      group=hoc_groups
      do while(group .gt. 0)
       point=hoc_points(group)
       do while(point .gt. 0)
        particle=hoc_particles(point)
        do while(particle .gt. 0)
         if(group .eq. 
     >     highest_level_group(particle_pointer(particle)))
     >     then
          part=particle_pointer(particle)
          d_x=pos_x(part)-x_0
          d_y=pos_y(part)-y_0
          d_z=pos_z(part)-z_0
          d2=d_x**2+d_y**2+d_z**2+
     >      4.0**(-level(group)-1)/float(grid_length**2)
          pot=-1.0/sqrt(d2)
          f_x=-d_x/d2**1.5
          f_y=-d_y/d2**1.5
          f_z=-d_z/d2**1.5
          write(98,98)group,level(group),particle,part,point,
     >      pos_x(part),pos_y(part),pos_z(part),
     >      f_x,force_x(part),
     >      f_y,force_y(part),f_z,force_z(part),
     >      pot,potential(part)
 98       format(5i8,11(1pe11.3))
         end if
         particle=next_particles(particle)
        end do
        point=next_points(point)
       end do
       group=next_groups(group)
      end do
c     
      return
      end

