      subroutine particle_lists(grid_length,level_max,level,hoc_groups,
     >  next_groups,point_pointer,pos_x,pos_y,pos_z,
     >  particle_pointer,hoc_particles,next_particles,
     >  pos_point_x,pos_point_y,pos_point_z,
     >  point_up_x,point_up_y,point_up_z,hoc_points,next_points,
     >  real_pointer,periodic,number_particles,
     >  groups_maxx,points_maxx,particles_maxx,particles_real)
c     
      implicit none
c     
      integer groups_maxx,points_maxx,particles_maxx,particles_real
c     
      integer length,grid_length,point,hoc_particles(points_maxx)
      integer p,number_particles,p_x,p_y,p_z,hoc_points(groups_maxx)
      integer next_points(points_maxx),next_groups(groups_maxx)
      integer particle_pointer(particles_maxx)
      integer next_particles(particles_maxx),new_particle,level_max
      integer lev,group,hoc_groups,particle
      integer pos_point_x(points_maxx),pos_point_y(points_maxx)
      integer pos_point_z(points_maxx),p_p
      integer point_up_x(points_maxx),point_up_y(points_maxx)
      integer point_up_z(points_maxx),level(groups_maxx)
      integer point_pointer(points_maxx),part,real_pointer(points_maxx)
c     
      real pos_x(particles_real),pos_y(particles_real)
      real pos_z(particles_real),grid_multiply,d_inv
      real a_grid_length
c     
      logical periodic
c     
      a_grid_length=grid_length
      grid_multiply=2**level_max*grid_length
c     
      if(periodic) then
        length=grid_length
      else
        length=grid_length+1
      end if
c     
      do point=1,length**3
        hoc_particles(point)=-1
      end do
c     
      do p=1,number_particles
        p_x=pos_x(p)*a_grid_length
        p_y=pos_y(p)*a_grid_length
        p_z=pos_z(p)*a_grid_length
        point=1+p_x+length*(p_y+length*p_z)
        particle_pointer(p)=p
        next_particles(p)=hoc_particles(point)
        hoc_particles(point)=p
      end do
c     
      new_particle=number_particles
c     
      if(level_max .lt. 1) return
c     
      do lev=1,level_max
        d_inv=1.0/float(2**(level_max-lev))
        group=hoc_groups
        do while(group .gt. 0)
          if(level(group) .eq. lev) then
            point=hoc_points(group)
            do while (point .gt. 0)
              hoc_particles(point)=-1
              point=next_points(point)
            end do
c     
            point=hoc_points(group)
            do while(point .gt. 0)
              p_p=point_pointer(point)
              if(real_pointer(point) .eq. 1) then
                particle=hoc_particles(p_p)
                do while(particle .gt. 0)
                  part=particle_pointer(particle)
                  p_x=(pos_x(part)*grid_multiply-
     >              pos_point_x(point))*d_inv
                  p_y=(pos_y(part)*grid_multiply-
     >              pos_point_y(point))*d_inv
                  p_z=(pos_z(part)*grid_multiply-
     >              pos_point_z(point))*d_inv
                  if(max(p_x,p_y,p_z) .gt. 1 .or.
     >              min(p_x,p_y,p_z) .lt. 0) stop 'die in lists'
                  p=point
                  if(p_x .gt. 0) p=point_up_x(p)
                  if(p_y .gt. 0) p=point_up_y(p)
                  if(p_z .gt. 0) p=point_up_z(p)
                  new_particle=new_particle+1
                  next_particles(new_particle)=hoc_particles(p)
                  hoc_particles(p)=new_particle
                  particle_pointer(new_particle)=part
                  particle=next_particles(particle)
                end do
              end if
              point=next_points(point)
            end do
          end if
          group=next_groups(group)
        end do
      end do
c     
 96   format(5i8,3(1pe14.5))
      return
      end
