      subroutine check_for_edge_trouble(pos_x,pos_y,pos_z,
     >  number_particles,particles_real)
c     
      implicit none
c     
      integer number_particles,particles_real,part
      real pos_x(particles_real),pos_y(particles_real)
      real pos_z(particles_real),epsilon
      data epsilon/5.0e-7/
c     
      do part=1,number_particles
       if(pos_x(part) .ge. 1.0) pos_x(part)=pos_x(part)-epsilon
       if(pos_y(part) .ge. 1.0) pos_y(part)=pos_y(part)-epsilon
       if(pos_z(part) .ge. 1.0) pos_z(part)=pos_z(part)-epsilon
c     
       if(pos_x(part) .le. 0.0) pos_x(part)=pos_x(part)+epsilon
       if(pos_y(part) .le. 0.0) pos_y(part)=pos_y(part)+epsilon
       if(pos_z(part) .le. 0.0) pos_z(part)=pos_z(part)+epsilon
      end do
c     
      return
      end
c
