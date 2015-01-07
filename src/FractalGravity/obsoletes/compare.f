      subroutine compare(group,hoc_groups,next_groups,
     >  hoc_points,next_points,point_pointer,level,level_max,
     >  pos_point_x,pos_point_y,pos_point_z,inside,
     >  potential_point,force_point_x,force_point_y,force_point_z,
     >  pos_x,pos_y,pos_z,potential,force_x,force_y,force_z,
     >  density,highest_level_group,number_particles,
     >  grid_length,particles_real,groups_maxx,points_maxx)
c     
      implicit none
c
      integer groups_maxx,points_maxx,particles_real
      integer group,hoc_groups,point,next_groups(groups_maxx)
      integer hoc_points(groups_maxx),level(10000)
      integer level_max,pos_point_x(points_maxx)
      integer grid_length,next_points(points_maxx)
      integer pos_point_y(points_maxx),pos_point_z(points_maxx)
      integer point_pointer(points_maxx),number_particles,part
      integer highest_level_group(particles_real)
c
      real potential_point(points_maxx),force_point_x(points_maxx)
      real force_point_y(points_maxx),force_point_z(points_maxx)
      real x_0,y_0,z_0,r_0,alpha,grid_multiply,eps,eps2,dx,dy,dz
      real cx,cy,cz,f_x,f_y,f_z,mass,pot,d,d2,x,y,z
      real potential(particles_real),force_x(particles_real)
      real force_y(particles_real),force_z(particles_real)
      real pos_x(particles_real),pos_y(particles_real)
      real pos_z(particles_real),density(points_maxx)
      real f_tot,force_tot,r
c
      logical inside(points_maxx)
c
      print*,'x_0,y_0,z_0,r_0,alpha'
      read(*,*)x_0,y_0,z_0,r_0,alpha
c
      grid_multiply=grid_length*2**level_max
      group=hoc_groups
      do while(group .gt. 0)
        point=hoc_points(group)
        eps=0.5/float(grid_length)/2.0**level(group)
        eps2=eps**2
        do while(point .gt. 0)
          x=pos_point_x(point)/grid_multiply
          y=pos_point_y(point)/grid_multiply
          z=pos_point_z(point)/grid_multiply
          dx=x-x_0
          dy=y-y_0
          dz=z-z_0
          d2=dx**2+dy**2+dz**2+eps2
          d=sqrt(d2)
          cx=dx/d
          cy=dy/d
          cz=dz/d
          if(d .gt. r_0) then
            f_x=-cx/d2
            f_y=-cy/d2
            f_z=-cz/d2
            pot=-1.0/d
          else
            if(alpha .gt. -1.99 .or. alpha .lt. -2.01) then
              mass=(d/r_0)**(3.0+alpha)
              f_x=-mass*cx/d2
              f_y=-mass*cy/d2
              f_z=-mass*cz/d2
              pot=-1.0/r_0-(1.0-(d/r_0)**(2.0+alpha))/(r_0*(2.0+alpha))
            else
              mass=d/r_0
              f_x=-mass*cx/d2
              f_y=-mass*cy/d2
              f_z=-mass*cz/d2
              pot=log(d/r_0)/r_0-1.0/r_0
            end if
          end if
c
          r=sqrt(dx**2+dy**2+dz**2)
          force_tot=sqrt(force_point_x(point)**2+
     >      force_point_y(point)**2+force_point_z(point)**2)
          f_tot=sqrt(f_x**2+f_y**2+f_z**2)
          write(99,1)level(group),group,inside(point),
     >      point,point_pointer(point),
     >      x,y,z,r,density(point),potential_point(point),pot,
     >      force_tot,f_tot,force_point_x(point),f_x,
     >      force_point_y(point),f_y,force_point_z(point),f_z
          point=next_points(point)
        end do
        group=next_groups(group)
      end do
c
 1    format(i2,i5,l2,i8,i8,15(1pe12.4))
c
      do part=1,number_particles
        group=highest_level_group(part)
        eps=0.5/float(grid_length)/2.0**level(group)
        eps2=eps**2
        dx=pos_x(part)-x_0
        dy=pos_y(part)-y_0
        dz=pos_z(part)-z_0
        d2=dx**2+dy**2+dz**2+eps2
        d=sqrt(d2)
        cx=dx/d
        cy=dy/d
        cz=dz/d
        if(d .gt. r_0) then
          f_x=-cx/d2
          f_y=-cy/d2
          f_z=-cz/d2
          pot=-1.0/d
        else
          if(alpha .gt. -1.99 .or. alpha .lt. -2.01) then
            mass=(d/r_0)**(3.0+alpha)
            f_x=-mass*cx/d2
            f_y=-mass*cy/d2
            f_z=-mass*cz/d2
            pot=-1.0/r_0-(1.0-(d/r_0)**(2.0+alpha))/(r_0*(2.0+alpha))
          else
            mass=d/r_0
            f_x=-mass*cx/d2
            f_y=-mass*cy/d2
            f_z=-mass*cz/d2
            pot=log(d/r_0)-1.0/r_0
          end if
        end if
        r=sqrt(dx**2+dy**2+dz**2)
        force_tot=sqrt(force_x(part)**2+
     >    force_y(part)**2+force_z(part)**2)
        f_tot=sqrt(f_x**2+f_y**2+f_z**2)
        write(98,2)part,group,level(group),
     >    pos_x(part),pos_y(part),pos_z(part),r,potential(part),pot,
     >    force_tot,f_tot,
     >    force_x(part),f_x,force_y(part),f_y,force_z(part),f_z
      end do
 2    format(i8,i5,i2,15(1pe12.4))
c
      return
      end
