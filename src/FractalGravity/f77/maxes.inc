      integer groups_max,particles_max,points_max
      integer particles_max_0,points_max_group
      integer high_points_max,double
      integer heavy,parts_heavy,particles_heavy
      parameter(double=1)
      parameter(groups_max=10000,particles_max=600000)
      parameter(particles_max_0=332000)
      parameter(points_max=700000,points_max_group=262143*double+1)
      parameter(high_points_max=points_max)
      parameter (heavy=0)
      parameter (parts_heavy=50000)
      parameter (particles_heavy=heavy*(parts_heavy-1)+1)
