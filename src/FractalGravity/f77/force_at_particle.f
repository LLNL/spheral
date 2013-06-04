      subroutine force_at_particle(group,
     >     hoc_points,next_points,potential_point,potential,
     >     force_point_x,force_point_y,force_point_z,particle_pointer,
     >     force_x,force_y,force_z,point_up_x,point_up_y,point_up_z,
     >     inside,hoc_particles,next_particles,pos_point_x,pos_point_y,
     >     pos_point_z,pos_x,pos_y,pos_z,level,level_max,grid_length,
     >     highest_level_group,debug,
     >     groups_maxx,points_maxx,particles_maxx,particles_real,tweaks)
c     
      implicit none
c     
      integer groups_maxx,points_maxx,particles_maxx,particles_real
      integer tweaks
      real force_point_x(points_maxx),force_point_y(points_maxx)
      real force_point_z(points_maxx)
      real force_x(particles_real),force_y(particles_real)
      real force_z(particles_real)
      real dens(8),f_x(8),f_y(8),f_z(8)
      real pott(8),d_inv
      real pos_x(particles_real),pos_y(particles_real)
      real pos_z(particles_real),scale,d_x,d_y,d_z,d_z_1
      real potential_point(points_maxx)
      real potential(particles_real),sum_prod,too_force_ful
c     
      logical inside(points_maxx),debug,bad_p,bad_x,bad_y,bad_z
      logical test_big
c     
      integer point,hoc_points(groups_maxx),particle
      integer level(groups_maxx)
      integer point_up_x(points_maxx),point_up_y(points_maxx)
      integer point_up_z(points_maxx)
      integer part,next_points(points_maxx)
      integer pos_point_x(points_maxx),pos_point_y(points_maxx)
      integer pos_point_z(points_maxx)
      integer p_x,p_y,p_z,p_x_y,p_z_x_y,p_z_x,p_z_y,grid_length
      integer level_max,hoc_particles(points_maxx)
      integer next_particles(particles_maxx),previous_point
      integer particle_pointer(particles_maxx)
      integer highest_level_group(particles_real),group
c     
      data too_force_ful/1.0e8/
c     
      previous_point=-10
      scale=grid_length*2**level_max
c     
c***  If a point has particles associated with it, then this point
c***  is at position (1,2,4,5,10,11,13,14)
c     
      d_inv=2.0**(level(group)-level_max)
      point=hoc_points(group)
      do while(point .gt. 0)
         particle=hoc_particles(point)
         do while(particle .gt. 0)
            part=particle_pointer(particle)
            if(group .eq.  highest_level_group(part)) then
               d_x=(pos_x(part)*scale-float(pos_point_x(point)))*d_inv
               d_y=(pos_y(part)*scale-float(pos_point_y(point)))*d_inv
               d_z=(pos_z(part)*scale-float(pos_point_z(point)))*d_inv
               if(abs(d_x-0.5) .gt. 0.5) stop 'die'
               if(abs(d_y-0.5) .gt. 0.5) stop 'die'
               if(abs(d_z-0.5) .gt. 0.5) stop 'die'
c     
               dens(1)=(1.0-d_x)*(1.0-d_y)
               dens(2)=d_x*(1.0-d_y)
               dens(3)=(1.0-d_x)*d_y
               dens(4)=d_x*d_y
c     
               d_z_1=1.0-d_z
c     
               dens(5)=dens(1)*d_z
               dens(6)=dens(2)*d_z
               dens(7)=dens(3)*d_z
               dens(8)=dens(4)*d_z
c     
               dens(1)=dens(1)*d_z_1
               dens(2)=dens(2)*d_z_1
               dens(3)=dens(3)*d_z_1
               dens(4)=dens(4)*d_z_1
c     
               if(previous_point .ne. point) then
                  p_x=point_up_x(point)
                  p_y=point_up_y(point)
                  p_x_y=point_up_y(p_x)
c     
                  pott(1)=potential_point(point)
                  f_x(1)=force_point_x(point)
                  f_y(1)=force_point_y(point)
                  f_z(1)=force_point_z(point)
c     
                  pott(2)=potential_point(p_x)
                  f_x(2)=force_point_x(p_x)
                  f_y(2)=force_point_y(p_x)
                  f_z(2)=force_point_z(p_x)
c     
                  pott(3)=potential_point(p_y)
                  f_x(3)=force_point_x(p_y)
                  f_y(3)=force_point_y(p_y)
                  f_z(3)=force_point_z(p_y)
c     
                  pott(4)=potential_point(p_x_y)
                  f_x(4)=force_point_x(p_x_y)
                  f_y(4)=force_point_y(p_x_y)
                  f_z(4)=force_point_z(p_x_y)
c     
                  p_z=point_up_z(point)
                  p_z_x=point_up_z(p_x)
                  p_z_y=point_up_z(p_y)
                  p_z_x_y=point_up_z(p_x_y)
c     
                  pott(5)=potential_point(p_z)
                  f_x(5)=force_point_x(p_z)
                  f_y(5)=force_point_y(p_z)
                  f_z(5)=force_point_z(p_z)
c     
                  pott(6)=potential_point(p_z_x)
                  f_x(6)=force_point_x(p_z_x)
                  f_y(6)=force_point_y(p_z_x)
                  f_z(6)=force_point_z(p_z_x)
c     
                  pott(7)=potential_point(p_z_y)
                  f_x(7)=force_point_x(p_z_y)
                  f_y(7)=force_point_y(p_z_y)
                  f_z(7)=force_point_z(p_z_y)
c     
                  pott(8)=potential_point(p_z_x_y)
                  f_x(8)=force_point_x(p_z_x_y)
                  f_y(8)=force_point_y(p_z_x_y)
                  f_z(8)=force_point_z(p_z_x_y)
c     
               end if
c     
               potential(part)=sum_prod(1,8,1,dens,pott)
               force_x(part)=sum_prod(1,8,1,dens,f_x)
               force_y(part)=sum_prod(1,8,1,dens,f_y)
               force_z(part)=sum_prod(1,8,1,dens,f_z)
c     
               bad_p=test_big(pott,1,8,1,too_force_ful)
               bad_x=test_big(f_x,1,8,1,too_force_ful)
               bad_y=test_big(f_y,1,8,1,too_force_ful)
               bad_z=test_big(f_z,1,8,1,too_force_ful)
               if(bad_p .or. bad_x .or. bad_y .or. bad_z) then
                  write(48,*)'bigg trouble ',group,point,particle,part
                  write(48,*)pos_x(part),pos_y(part),pos_z(part)
                  write(48,*)pos_point_x(point),pos_point_y(point),
     >                 pos_point_z(point)
                  write(48,*)pott
                  write(48,*)f_x
                  write(48,*)f_y
                  write(48,*)f_z
               end if
            end if
c     
            particle=next_particles(particle)
         end do
c     
         previous_point=point
         point=next_points(point)
      end do
      return
      end
c
