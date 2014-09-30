      subroutine negatives_particles(particles_maxx,particles_real,
     >  next_particles,
     >  particle_pointer,highest_level_group,
     >  potential,force_x,force_y,force_z,tweaks)
c     
      implicit none
c     
      integer particles_maxx,particle,particles_real,tweaks
      integer next_particles(particles_maxx)
      integer particle_pointer(particles_maxx)
      integer highest_level_group(particles_real)
c     
      real potential(particles_real),force_x(particles_real)
      real force_y(particles_real),force_z(particles_real)

c     
      do particle=1,particles_maxx
       next_particles(particle)=-particle
       particle_pointer(particle) =-particle
      end do
c     
      do particle=1,particles_real
       force_x(particle)=-1.0e30*float(particle)
       force_y(particle)=-1.0e30*float(particle)
       force_z(particle)=-1.0e30*float(particle)
       potential(particle)=-1.0e30*float(particle)
       highest_level_group(particle)=-particle
      end do
c     
      return
      end
