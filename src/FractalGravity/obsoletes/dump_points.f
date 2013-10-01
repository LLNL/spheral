      subroutine dump_points(hoc_groups,next_groups,level,
     >  level_max,hoc_points,next_points,
     >  pos_point_x,pos_point_y,pos_point_z,
     >  pos_x,pos_y,pos_z,particle_pointer,inside,
     >  highest_level_group,hoc_particles,next_particles,
     >  number_particles,groups_maxx,points_maxx,particles_maxx)
c     
      implicit none
c     
      integer groups_maxx,points_maxx,particles_maxx,number_particles
      integer hoc_groups,next_groups(groups_maxx),level(groups_maxx)
      integer level_max,hoc_points(groups_maxx),next_points(points_maxx)
      integer pos_point_x(points_maxx),pos_point_y(points_maxx)
      integer pos_point_z(points_maxx)
      integer highest_level_group(particles_maxx)
      integer particle_pointer(particles_maxx)
      integer next_particles(particles_maxx),hoc_particles(points_maxx)
      integer lev,group,point,part
      real pos_x(particles_maxx),pos_y(particles_maxx)
      real pos_z(particles_maxx)
c
      logical inside(points_maxx)
c     
      do lev=0,level_max
        open(unit=lev+10,form='formatted')
        open(unit=lev+20,form='formatted')
      end do
c
      group=hoc_groups
      do while(group .gt. 0)
        lev=level(group)
        point=hoc_points(group)
        do while(point .gt. 0)
          write(10+lev,10)group,point,inside(point),
     >      pos_point_x(point),pos_point_y(point),pos_point_z(point)
          point=next_points(point)
        end do
        group=next_groups(group)
      end do
 10   format(i5,i8,l2,3i8)
c
      do part=1,number_particles
        group=highest_level_group(part)
        lev=level(group)
        write(20+lev,20)group,part,
     >    pos_x(part),pos_y(part),pos_z(part)
      end do
c
 20   format(i5,i8,3(1pe13.5))
      return
      end
