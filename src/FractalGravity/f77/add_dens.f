      subroutine add_dens(dens,dm,d_x,d_y,d_z)
c
      implicit none
c
      real dens(8),d_x,d_y,d_z,d_1,d_2,d_3,d_4,d_z_1,dm
c
      d_1=(1.0-d_x)*(1.0-d_y)
      d_2=d_x*(1.0-d_y)
      d_3=(1.0-d_x)*d_y
      d_4=d_x*d_y
      d_z_1=1.0-d_z
c
      dens(1)=dens(1)+d_1*d_z_1*dm
      dens(2)=dens(2)+d_2*d_z_1*dm
      dens(3)=dens(3)+d_3*d_z_1*dm
      dens(4)=dens(4)+d_4*d_z_1*dm
      dens(5)=dens(5)+d_1*d_z*dm
      dens(6)=dens(6)+d_2*d_z*dm
      dens(7)=dens(7)+d_3*d_z*dm
      dens(8)=dens(8)+d_4*d_z*dm
c
      return
      end
