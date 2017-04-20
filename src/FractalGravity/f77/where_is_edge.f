      subroutine where_is_edge(highest_level_group,
     >     mother_group,highest_level_point,number_particles,
     >     fake_particle,point_pointer,
     >     real_pointer,up_x,up_y,up_z,down_x,down_y,down_z,ins,
     >     number_fakes,groupies,long,long_edge,
     >     particles_real)
c     
      implicit none
c     
      integer particles_real,groupies
      integer long,long_edge,number_fakes
      integer point,real_pointer(long),up_x(long),up_y(long),up_z(long)
      integer down_x(long),down_y(long),down_z(long),p,p0
      integer rp,p1,p2,p3,particle,point_pointer(long)
      integer highest_level_point(particles_real),number_particles
      integer fake_particle(long_edge),highest_level_group(long)
      integer mother_group(groupies)
c     
      logical ins(long),really_inside
c     
      number_fakes=0
c     
      do particle=1,number_particles
         if(highest_level_group(particle) .gt. 1) then
            point=highest_level_point(particle)
            rp=real_pointer(point)-1
            p3=rp/9
            p2=mod(rp/3,3)
            p1=mod(rp,3)
            p=point
            if(p1 .eq. 1) p=up_x(p)
            if(p2 .eq. 1) p=up_y(p)
            if(p3 .eq. 1) p=up_z(p)
            p0=p
c     
            if(.not. ins(p)) then
               p=point
               if(p1 .eq. 1) p=down_x(p)
               if(p2 .eq. 1) p=down_y(p)
               if(p3 .eq. 1) p=down_z(p)
               highest_level_point(particle)=point_pointer(p)
               highest_level_group(particle)=
     >              mother_group(highest_level_group(particle))
c               write(91,*)'on the edge ',particle,point,p,p0,
c     >              real_pointer(point),real_pointer(p),real_pointer(p0)
            else
c     
               p=point
               if(p1 .eq. 1) p=down_x(p)
               if(p2 .eq. 1) p=down_y(p)
               if(p3 .eq. 1) p=down_z(p)
c     
               really_inside=.false.
               if(ins(p)) then
                  p=up_x(up_x(p))
                  if(ins(p)) then
                     p=up_y(up_y(p))
                     if(ins(p)) then
                        p=down_x(down_x(p))
                        if(ins(p)) then
                           p=up_z(up_z(p))
                           if(ins(p)) then
                              p=up_x(up_x(p))
                              if(ins(p)) then
                                 p=down_y(down_y(p))
                                 if(ins(p)) then
                                    p=down_x(down_x(p))
                                    if(ins(p)) then
                                       really_inside=.true.
                                    end if
                                 end if
                              end if
                           end if
                        end if
                     end if
                  end if
               end if
               if(.not. really_inside) then
                  number_fakes=number_fakes+1
                  fake_particle(number_fakes)=particle
c               write(92,*)'not on the edge ',particle,point,p,
c     >              real_pointer(point),real_pointer(p)
               end if
            end if
         end if
      end do
      return
      end
