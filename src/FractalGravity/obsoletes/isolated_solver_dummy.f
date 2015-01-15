      subroutine isolated_solver(grid_length,density,
     >  potential_point,debug,memory_value,
     >  groups_maxx,points_maxx,particles_maxx,tweaks)
c
      implicit none
c
      integer grid_length,groups_maxx,points_maxx
      integer particles_maxx,tweaks,memory_value
      logical debug
      real density(*),potential_point(*)
c
      stop 'isolated dummy is never supposed to be called'
      return
      end
