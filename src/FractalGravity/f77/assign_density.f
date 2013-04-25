      subroutine assign_density(group,hoc_particles,density,
     >     next_particles,grid_length,level_max,
     >     pos_x,pos_y,pos_z,particle_mass,density_0,
     >     pos_point_x,pos_point_y,pos_point_z,level,
     >     number_particles,periodic,hoc_points,next_points,
     >     particle_pointer,point_up_x,point_up_y,point_up_z,
     >     point_down_x,point_down_y,point_down_z,debug,
     >     point_pointer,inside,set_zero,set_dens,set_scaling,
     >     groups_maxx,points_maxx,particles_maxx,
     >     particles_real,tweaks)
c     
      implicit none
c     
      integer groups_maxx,points_maxx,particles_maxx,particles_real
      integer hoc_particles(points_maxx),tweaks
      integer point,point_pointer(points_maxx)
      integer next_particles(particles_maxx),grid_length,level_max
      integer pos_point_x(points_maxx),pos_point_y(points_maxx)
      integer pos_point_z(points_maxx),part,particle
      integer level(groups_maxx),hoc_points(points_maxx)
      integer next_points(points_maxx)
      integer particle_pointer(particles_maxx)
      integer point_up_x(points_maxx),point_up_y(points_maxx)
      integer point_up_z(points_maxx)
      integer group,number_particles
      integer point_down_x(points_maxx),point_down_y(points_maxx)
      integer point_down_z(points_maxx),p_p
c     
      integer p_up_x,p_up_y,p_up_x_y,p_up_z,p_down_x,p_down_y,p_down_z
      logical periodic,debug,inside(points_maxx)
      logical set_dens,set_scaling,set_zero
c     
      real dens(8),pos_x(particles_real),pos_y(particles_real)
      real pos_z(particles_real),particle_mass(particles_real)
      real scale,d_x,d_y,d_z,density(points_maxx),d_inv
      real scaling,density_0
c     
      scale=grid_length*2**level_max
      d_inv=2.0**(level(group)-level_max)
c     
      if(set_zero) then
         point=hoc_points(group)
         do while(point .gt. 0)
            density(point)=0.0
            point=next_points(point)
         end do
      end if
c     
      point=hoc_points(group)
      do while(point .gt. 0)
c***  particles_at_point(point) is positive only at the lower left point 
c***  in a virtual cell
c***  the density at a point that is not inside may be incorrect.
c***  that is OK, we only use the density at inside points.
         particle=hoc_particles(point)
         if(particle .gt. 0) then
            call set_constant_real(dens,1,8,1,0.0)
            do while(particle .gt. 0)
               part=particle_pointer(particle)
               if(particle_mass(part) .ne. 0.0) then
c     
                  d_x=(pos_x(part)*scale-float(pos_point_x(point)))*
     >                 d_inv
                  d_y=(pos_y(part)*scale-float(pos_point_y(point)))*
     >                 d_inv
                  d_z=(pos_z(part)*scale-float(pos_point_z(point)))*
     >                 d_inv
                  if(
     >                 abs(d_x-0.5) .gt. 0.5 .or.
     >                 abs(d_y-0.5) .gt. 0.5 .or.
     >                 abs(d_z-0.5) .gt. 0.5) then
                     write(44,*)'die ',group,point,part,level(group),
     >                    d_inv,
     >                    d_x,d_y,d_z,
     >                    pos_x(part),pos_y(part),pos_z(part),
     >                    pos_point_x(point),
     >                    pos_point_y(point),pos_point_z(point)
                  end if
c     
                  call add_dens(dens,particle_mass(part),d_x,d_y,d_z)
c     
               end if
               particle=next_particles(particle)
            end do
c     
            p_up_x=point_up_x(point)
            p_up_y=point_up_y(point)
            p_up_x_y=point_up_x(p_up_y)
c     
            density(point)=density(point)+dens(1)
            density(p_up_x)=density(p_up_x)+dens(2)
            density(p_up_y)=density(p_up_y)+dens(3)
            density(p_up_x_y)=density(p_up_x_y)+dens(4)
c     
            density(point_up_z(point))=
     >           density(point_up_z(point))+dens(5)
            density(point_up_z(p_up_x))=
     >           density(point_up_z(p_up_x))+dens(6)
            density(point_up_z(p_up_y))=
     >           density(point_up_z(p_up_y))+dens(7)
            density(point_up_z(p_up_x_y))=
     >           density(point_up_z(p_up_x_y))+dens(8)
c***  does this point always exist???? It must!!
         end if
         point=next_points(point)
      end do
c     
      if(set_scaling) then
         scaling=(float(grid_length)*2.0**level(group))**3
         point=hoc_points(group)
         do while(point .gt. 0)
            if(inside(point)) then
               density(point)=density(point)*scaling
            else if(point_pointer(point) .gt. 0) then
               density(point)=density(point_pointer(point))
            end if
            point=next_points(point)
         end do
c     
         if(group .gt. 1) then
            point=hoc_points(group)
            do while(point .gt. 0)
               if((.not. inside(point)) .and. 
     >              point_pointer(point) .lt. 0) then
                  p_up_x=point_up_x(point)
                  p_up_y=point_up_y(point)
                  p_up_z=point_up_z(point)
                  p_down_x=point_down_x(point)
                  p_down_y=point_down_y(point)
                  p_down_z=point_down_z(point)
                  p_p=-point_pointer(point)
c     
                  if(debug)
     >                 write(44,*)group,point,p_p,p_up_x,p_up_y,p_up_z,
     >                 p_down_x,p_down_y,p_down_z
c     
                  if(p_p .eq. 2 .or. p_p .eq. 8 .or.
     >                 p_p .eq. 20 .or. p_p .eq. 26) then
                     density(point)=0.5*(
     >                    density(p_up_x)+
     >                    density(p_down_x))
c     
                  else if(p_p .eq. 4 .or. p_p .eq. 6 .or. 
     >                    p_p .eq. 22 .or. p_p .eq. 24) then
                     density(point)=0.5*(density(p_up_y)+
     >                    density(p_down_y))
c     
                  else if(p_p .eq. 10 .or. p_p .eq. 12 .or. 
     >                    p_p .eq. 16 .or. p_p .eq. 18) then
                     density(point)=0.5*(density(p_up_z)+
     >                    density(p_down_z))
c     
                  else if(p_p .eq. 5 .or. p_p .eq. 23) then
                     density(point)=0.25*(
     >                    density(point_down_x(p_down_y))+
     >                    density(point_up_x(p_down_y))+
     >                    density(point_down_x(p_up_y))+
     >                    density(point_up_x(p_up_y)))
                  else if(p_p .eq. 11 .or. p_p .eq. 17) then
                     density(point)=0.25*(
     >                    density(point_down_x(p_down_z))+
     >                    density(point_up_x(p_down_z))+
     >                    density(point_down_x(p_up_z))+
     >                    density(point_up_x(p_up_z)))
                  else if(p_p .eq. 13 .or. p_p .eq. 15) then
                     density(point)=0.25*(
     >                    density(point_down_y(p_down_z))+
     >                    density(point_up_y(p_down_z))+
     >                    density(point_down_y(p_up_z))+
     >                    density(point_up_y(p_up_z)))
                  else if(p_p .eq. 14) then
                     print*,p_p
                     stop 'density not with you at point '
                  else if(p_p .ge. 0) then
                     print*,p_p
                     stop 'density not with you at point '
                  end if
               end if
               point=next_points(point)
            end do
         end if
      end if
c     
      if(set_dens) then
         point=hoc_points(group)
         do while(point .gt. 0)
            if(inside(point)) density(point)=density(point)-density_0
            point=next_points(point)
         end do
      end if
      return
      end
