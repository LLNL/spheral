      subroutine tree_start(grid_length,level_max,group,groups_maxx,
     >  points_maxx,particles_maxx,periodic,hoc_points,
     >  next_points,highest_used_point,highest_used_particle,
     >  point_up_x,point_up_y,point_up_z,point_pointer,real_pointer,
     >  hoc_particles,next_particles,particle_pointer,
     >  point_down_x,point_down_y,point_down_z,pos_point_x,
     >  pos_point_y,pos_point_z,inside,particles_at_point,moat,
     >  number_particles,particles_real,pos_x,pos_y,pos_z,debug,
     >  tweaks)
c     
      implicit none
c     
      integer number_particles,particles_real,particle,tweaks
      integer n,i,j,k,g,grid_length,length,length_2,length_3,n_p
      integer grid_x,grid_y,grid_z,grid_multiply,level_max
      integer points_maxx,groups_maxx,group,particles_maxx,point
      integer x_up,y_up,z_up,x_down,y_down,z_down,memory_value
      integer point_x,point_y,point_z,moat,highest_used_particle
      integer hoc_points(groups_maxx),next_points(points_maxx)
      integer point_up_x(points_maxx),point_down_x(points_maxx)
      integer point_up_y(points_maxx),point_down_y(points_maxx)
      integer point_up_z(points_maxx),point_down_z(points_maxx)
      integer pos_point_x(points_maxx),pos_point_y(points_maxx)
      integer pos_point_z(points_maxx),particles_at_point(points_maxx)
      integer highest_used_point,particle_pointer(particles_maxx)
      integer point_pointer(points_maxx),real_pointer(points_maxx)
      integer hoc_particles(points_maxx),next_particles(particles_maxx)
c     
      real pos_x(particles_real),pos_y(particles_real)
      real pos_z(particles_real),a_grid_length
c     
      logical periodic,debug,inside(points_maxx)
c     
      n(i,j,k,g)=i+(j-1+(k-1)*g)*g
c     
c***  call CHECK_FOR_EDGE_TROUBLE
c***  checking for round off errors causing trouble for wraparound
c
      if(debug)print*,'treestart'
c     
      length=grid_length
      if(.not. periodic) length=grid_length+1
      a_grid_length=grid_length
      grid_multiply=2**level_max
      length_2=length**2
      length_3=length**3
c     
      memory_value=0
      if(length_3 .ge. points_maxx) then
        memory_value=1
        print*,'tree_start memory value= ',memory_value,length_3,
     >    points_maxx
        return
      end if
c     
      hoc_points(group)=-1
      do grid_z=1,length
        do grid_y=1,length
          do grid_x=1,length
            point=n(grid_x,grid_y,grid_z,length)
            point_pointer(point)=point
            real_pointer(point)=1
            next_points(point)=hoc_points(group)
            hoc_points(group)=point
            inside(point)=.true.
            hoc_particles(point)=-1
            particles_at_point(point)=0
            pos_point_x(point)=(grid_x-1)*grid_multiply
            pos_point_y(point)=(grid_y-1)*grid_multiply
            pos_point_z(point)=(grid_z-1)*grid_multiply
          end do
        end do
      end do
c     
      highest_used_point=n(length,length,length,length)
c     
      if(periodic) then
        do grid_z=1,length
          do grid_y=1,length
            do grid_x=1,length
              x_up=mod(grid_x,length)+1
              y_up=mod(grid_y,length)+1
              z_up=mod(grid_z,length)+1
              x_down=mod(grid_x+length-2,length)+1
              y_down=mod(grid_y+length-2,length)+1
              z_down=mod(grid_z+length-2,length)+1
              point=n(grid_x,grid_y,grid_z,length)
              point_up_x(point)=n(x_up,grid_y,grid_z,length)
              point_up_y(point)=n(grid_x,y_up,grid_z,length)
              point_up_z(point)=n(grid_x,grid_y,z_up,length)
              point_down_x(point)=n(x_down,grid_y,grid_z,length)
              point_down_y(point)=n(grid_x,y_down,grid_z,length)
              point_down_z(point)=n(grid_x,grid_y,z_down,length)
            end do
          end do
        end do
      else
        do grid_z=1,length
          do grid_y=1,length
            do grid_x=1,length
              point=n(grid_x,grid_y,grid_z,length)
              x_up=grid_x+1
              y_up=grid_y+1
              z_up=grid_z+1
              if(x_up .gt. length) then
                point_up_x(point)=-1
                inside(point)=.false.
              else
                point_up_x(point)=n(x_up,grid_y,grid_z,length)
              end if
              if(y_up .gt. length) then
                point_up_y(point)=-1
                inside(point)=.false.
              else
                point_up_y(point)=n(grid_x,y_up,grid_z,length)
              end if
              if(z_up .gt. length) then
                point_up_z(point)=-1
                inside(point)=.false.
              else
                point_up_z(point)=n(grid_x,grid_y,z_up,length)
              end if
c     
              x_down=grid_x-1
              y_down=grid_y-1
              z_down=grid_z-1
              if(x_down .lt. 1) then
                point_down_x(point)=-1
                inside(point)=.false.
              else
                point_down_x(point)=n(x_down,grid_y,grid_z,length)
              end if
              if(y_down .lt. 1) then
                point_down_y(point)=-1
                inside(point)=.false.
              else
                point_down_y(point)=n(grid_x,y_down,grid_z,length)
              end if
              if(z_down .lt. 1) then
                point_down_z(point)=-1
                inside(point)=.false.
              else
                point_down_z(point)=n(grid_x,grid_y,z_down,length)
              end if
            end do
          end do
        end do
      end if
c     
      if(periodic) then
        do particle=1,number_particles
          if(pos_x(particle) .ge. 1.0) then
            n_p=pos_x(particle)
            pos_x(particle)=pos_x(particle)-float(n_p)
          else if(pos_x(particle) .lt. 0.0) then
            n_p=-pos_x(particle)
            pos_x(particle)=pos_x(particle)+float(n_p+1)
          end if
c     
          if(pos_y(particle) .ge. 1.0) then
            n_p=pos_y(particle)
            pos_y(particle)=pos_y(particle)-float(n_p)
          else if(pos_y(particle) .lt. 0.0) then
            n_p=-pos_y(particle)
            pos_y(particle)=pos_y(particle)+float(n_p+1)
          end if
c     
          if(pos_z(particle) .ge. 1.0) then
            n_p=pos_z(particle)
            pos_z(particle)=pos_z(particle)-float(n_p)
          else if(pos_z(particle) .lt. 0.0) then
            n_p=-pos_z(particle)
            pos_z(particle)=pos_z(particle)+float(n_p+1)
          end if
        end do
c     
        call check_for_edge_trouble(pos_x,pos_y,pos_z,
     >    number_particles,particles_real)
c     
      end if
c     
      highest_used_particle=0
c     
      do particle=1,number_particles
        point_x=pos_x(particle)*a_grid_length
        point_y=pos_y(particle)*a_grid_length
        point_z=pos_z(particle)*a_grid_length
        point=n(point_x+1,point_y+1,point_z+1,length)
        if(periodic) then
          next_particles(particle)=hoc_particles(point)
          hoc_particles(point)=particle
          particles_at_point(point)=particles_at_point(point)+1
          particle_pointer(particle)=particle
        else
          if(point_x .ge. moat .and. point_y .ge. moat .and. 
     >      point_z .ge. moat .and.
     >      point_x .lt. grid_length-moat .and.
     >      point_y .lt. grid_length-moat .and.
     >      point_z .lt. grid_length-moat) then
            next_particles(particle)=hoc_particles(point)
            hoc_particles(point)=particle
            particles_at_point(point)=particles_at_point(point)+1
            particle_pointer(particle)=particle
          else
            particle_pointer(particle)=-particle
          end if              
        end if
      end do
c     
      highest_used_particle=number_particles
c     
c***  The moat is incorporated so that at level=0 no particle
c***  is in a virtual cell that borders on vacuum. This is only a problem
c***  with isolated boundary conditions
c     
c     
      return
      end
