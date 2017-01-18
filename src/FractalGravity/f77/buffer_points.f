      subroutine buffer_points(group,really_high,n_highs,
     >   it_is_high,level,grid_length,level_max,
     >   pos_x,pos_y,pos_z,hoc_points,next_points,
     >   pos_point_x,pos_point_y,pos_point_z,
     >   point_up_x,point_up_y,point_up_z,
     >   point_down_x,point_down_y,point_down_z,hoc_particles,
     >   next_particles,particle_pointer,minimum_number,
     >   periodic,inside,padding,
     >   groups_maxx,points_maxx,particles_maxx,particles_real)
c     
      implicit none
c     
      include 'maxx.inc'
c     
      integer points_maxx,groups_maxx,particles_maxx,particles_real
      integer high_points_maxx_group
c     
      integer group,point,corner,numbers(8),hoc_points(groups_maxx)
      integer n,n_highs,really_high(high_points_max_group)
      integer point_up_x(points_maxx),point_up_y(points_maxx)
      integer point_up_z(points_maxx),point_down_x(points_maxx)
      integer point_down_y(points_maxx),point_down_z(points_maxx)
      integer pos_point_x(points_maxx),pos_point_y(points_maxx)
      integer pos_point_z(points_maxx),particle_pointer(particles_maxx)
      integer grid_multiply,grid_length,level_max,level(groups_maxx)
      integer n_x,n_y,n_z,particle,part,minimum_number,padding
      integer hoc_particles(points_maxx),next_particles(particles_maxx)
      integer next_points(points_maxx),point_x,point_y,point_z
      integer pad_x,pad_y,pad_z,p,extras(7)
c     
      real pos_x(particles_real),pos_y(particles_real)
      real pos_z(particles_real),d_inv
c     
      logical it_is_high(-1:points_maxx),check_high,periodic
      logical inside(points_maxx),careful
c     
c***  call LIST_BUFFER
c***  add the buffer points
c     
      careful=.not. periodic .and. group .eq. 1 .and. padding .gt. 1
      if(careful) stop 'big padding not yet implemented'
c     
      high_points_maxx_group=high_points_max_group
      grid_multiply=grid_length*2**level_max
      d_inv=2.0**(-(level_max-level(group)-1))
c
      if(padding .gt. 0) then
        do n=1,n_highs
          point=really_high(n)
          do p=1,padding
            point=point_down_x(point)
          end do
          do p=1,padding
            point=point_down_y(point)
          end do
          do p=1,padding
            point=point_down_z(point)
          end do
c     
          point_z=point
          do pad_z=-padding,padding
            point_y=point_z
            do pad_y=-padding,padding
              point_x=point_y
              do pad_x=-padding,padding
                 it_is_high(point_x)=.true.
                point_x=point_up_x(point_x)
              end do
              point_y=point_up_y(point_y)
            end do
            point_z=point_up_z(point_z)
          end do
        end do
      else
        do n=1,n_highs
          point=really_high(n)
c     
          do corner=1,8
            numbers(corner)=0
          end do
c     
          particle=hoc_particles(point)
          do while(particle .gt. 0)
            part=particle_pointer(particle)
            n_x=(pos_x(part)*grid_multiply-pos_point_x(point))*d_inv
            n_y=(pos_y(part)*grid_multiply-pos_point_y(point))*d_inv
            n_z=(pos_z(part)*grid_multiply-pos_point_z(point))*d_inv
c     
            if(n_x .lt. 0 .or. n_x .gt. 1) stop 'die in buffer points x'
            if(n_y .lt. 0 .or. n_y .gt. 1) stop 'die in buffer points y'
            if(n_z .lt. 0 .or. n_z .gt. 1) stop 'die in buffer points z'
c     
            corner=1+n_x+2*n_y+4*n_z
            numbers(corner)=numbers(corner)+1
c     
            particle=next_particles(particle)
          end do
c
          do corner=1,8
            if(check_high(corner,numbers,minimum_number,8)) then
              call list_buffer(point,point_up_x,point_up_y,point_up_z,
     >           point_down_x,point_down_y,point_down_z,corner,
     >           it_is_high,extras,points_maxx)
              do p=1,7
                 it_is_high(extras(p))=.true.
              end do
            end if
          end do
        end do
      end if
c     
      if(group .eq. 1 .and. .not. periodic) then
        point=hoc_points(group)
        do while(point .gt. 0)
          it_is_high(point)=it_is_high(point) .and. inside(point)
          point=next_points(point)
        end do
      end if
c     
      return
      end
