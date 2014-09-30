      subroutine find_highest_level_group(
     >     hoc_groups,next_groups,
     >     hoc_points,next_points,hoc_particles,
     >     next_particles,particle_pointer,number_particles,
     >     heavy_version,force_max,mother_group,
     >     particle_mass,grid_length,down_x,down_y,down_z,
     >     real_pointer,point_pointer,highest_level_point,
     >     level,highest_level_group,debug,
     >     groups_maxx,points_maxx,particles_maxx,
     >     particles_real,tweaks)
c     
      implicit none
c     
      logical heavy_version,debug
c     
      integer groups_maxx,points_maxx,particles_maxx
      integer particles_real,tweaks,h_p,grid_length
      integer group,point,hoc_points(groups_maxx),particle,pp
      integer next_points(points_maxx),hoc_particles(points_maxx)
      integer highest_level_group(particles_real)
      integer next_particles(particles_maxx),part
      integer particle_0,hoc_groups,next_groups(groups_maxx)
      integer particle_pointer(particles_maxx),level(groups_maxx)
      integer p,number_particles,mother_group(groups_maxx)
      integer real_pointer(points_maxx),point_pointer(points_maxx)
      integer highest_level_point(particles_real)
      integer down_x(particles_maxx),down_y(particles_maxx)
      integer down_z(particles_maxx)
c     
      real force_max,log4,tmp_mass,very_small
      real particle_mass(particles_real)
c     
      if(debug) print*,'enter find highest level group'
c     
      do p=1,number_particles
         highest_level_group(p)=1
      end do
c     
      point=hoc_points(1)
      do while(point .gt. 0)
         particle=hoc_particles(point)
         do while(particle .gt. 0)
            highest_level_point(particle)=point
            particle=next_particles(particle)
         end do
         point=next_points(point)
      end do
c     
      group=hoc_groups
      do while(group .gt. 0)
         if(group .gt. 1) then
            point=hoc_points(group)
            do while(point .gt. 0)
               particle=hoc_particles(point)
               do while(particle .gt. 0)
                  particle_0=particle_pointer(particle)
                  if(level(group).gt.
     >                 level(highest_level_group(particle_0))) then
                     highest_level_group(particle_0)=group
                     highest_level_point(particle_0)=point
                  else if(level(group).eq.
     >                    level(highest_level_group(particle_0))) then
                     write(32,*)particle,particle_0,group,
     >                    highest_level_group(particle_0),level(group),
     >                    level(highest_level_group(particle_0)),point
                     stop 'particle in 2 groups at same level'
                  end if
                  particle=next_particles(particle)
               end do
               point=next_points(point)
            end do
         end if
         group=next_groups(group)
      end do
c     
      if(heavy_version) then
         log4=alog(4.0)
         tmp_mass=force_max/float(grid_length**2)
         very_small=1.0e-30
         do part=1,number_particles
            h_p=alog(tmp_mass/(particle_mass(part)+very_small))/log4
            do while(level(highest_level_group(part)) .gt. h_p)
               highest_level_group(part)=
     >              mother_group(highest_level_group(part))
               p=highest_level_point(p)
               pp=real_pointer(p)-1
               if(pp/9 .eq. 1) p=down_z(p)
               if(mod(pp/3,3) .eq. 1) p=down_y(p)
               if(mod(pp,3) .eq. 1) p=down_x(p)
               highest_level_point(part)=point_pointer(p)
            end do
         end do
      end if

      if(debug) print*,'exit find highest level group'
c     
      return
      end

