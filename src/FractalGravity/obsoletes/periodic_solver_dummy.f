      subroutine periodic_solver(grid_length,number_particles,
     >  density,potential_point,debug,
     >  groups_maxx,points_maxx,particles_maxx,tweaks)
c     
      implicit none
      integer grid_length,number_particles
      integer groups_maxx,points_maxx,particles_maxx,tweaks
      real density,potential_point
      logical debug
c
      stop 'periodic dummy is never supposed to be called'
      return
      end
