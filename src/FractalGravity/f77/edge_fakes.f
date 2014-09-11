      subroutine edge_fakes(grid_length,pos_x,pos_y,pos_z,
     >     pos_x_tmp,pos_y_tmp,pos_z_tmp,
     >     particle_mass,particle_mass_tmp,
     >     level,number_particles,number_particles_tmp,
     >     highest_level_point,point_down_x,point_down_y,
     >     point_down_z,it_is_high,particles_at_point,
     >     point_pointer,real_pointer,minimum_number,
     >     highest_level_group,level_min,level_max,
     >     particles_heavy_tmp,groupies,long,particles_real)
c     
      implicit none
c     
      integer particles_real,groupies,grid_length,long
      integer lev,level_max,p,number_particles,h_p,level_min
      integer highest_level_group(particles_real),p_tmp,n_p,p_x,p_y,p_z
      integer level(groupies),number_particles_tmp,pp
      integer particles_heavy_tmp,rp,rp1,rp2,rp3,point
      integer highest_level_point(long),point_down_x(long)
      integer point_down_y(long)
      integer point_down_z(long),particles_at_point(long)
      integer real_pointer(long),minimum_number
      integer point_pointer(long)
c     
      logical check_high,faker,it_is_high(long)
c     
      real spacing(0:8)
      real delta,off,small,offset(0:8),smaller(0:8)
      real pos_x(particles_real),pos_y(particles_real)
      real pos_z(particles_real),particle_mass(particles_real)
      real pos_x_tmp(particles_heavy_tmp),pos_y_tmp(particles_heavy_tmp)
      real pos_z_tmp(particles_heavy_tmp)
      real particle_mass_tmp(particles_heavy_tmp)
c     
      do lev=0,level_max
         spacing(lev)=1.0/float(grid_length*2**lev)
         offset(lev)=-0.5*spacing(lev)
         smaller(lev)=8.0**(-lev)
      end do
c     
      p_tmp=number_particles_tmp
      do p=1,number_particles
c     if(particle_mass(p) .gt. 1.0e-7 .and. 
c     >        highest_level_group(p) .gt. 1) then
c     point=highest_level_point(p)
c     write(96,*)p,point,particles_at_point(point),
c     >           highest_level_group(p),it_is_high(point),
c     >           pos_x(p),pos_y(p),pos_z(p)
c     end if
         lev=level(highest_level_group(p))
         if(lev .gt. 0) then
            point=highest_level_point(p)
            pp=point
            rp=real_pointer(pp)-1
            rp3=rp/9
            rp2=mod(rp/3,3)
            rp1=mod(rp,3)
            if(rp3 .eq. 1) pp=point_down_z(pp)
            if(rp2 .eq. 1) pp=point_down_y(pp)
            if(rp1 .eq. 1) pp=point_down_x(pp)
            pp=point_pointer(pp)
            faker=it_is_high(pp) .and. .not.
     >           check_high(pp,particles_at_point,
     >           minimum_number,long)
            if(faker) then
               h_p=lev-1
               do while(faker .and. h_p .gt. 0)
                  h_p=h_p-1
                  rp=real_pointer(pp)-1
                  rp3=rp/9
                  rp2=mod(rp/3,3)
                  rp1=mod(rp,3)
                  if(rp3 .eq. 1) pp=point_down_z(pp)
                  if(rp2 .eq. 1) pp=point_down_y(pp)
                  if(rp1 .eq. 1) pp=point_down_x(pp)
                  pp=point_pointer(pp)
                  faker=it_is_high(pp) .and. h_p .gt. 0 .and. .not.
     >                 check_high(pp,particles_at_point,
     >                 minimum_number,long)
               end do
c     
               n_p=2**(lev-h_p)
               if(p_tmp+n_p**3+1 .gt. particles_heavy_tmp) 
     >              stop 'too many fakers'
               delta=spacing(h_p)/float(n_p)
               off=-delta*float(n_p-1)*0.5
               small=smaller(lev-h_p)
               do p_z=1,n_p
                  do p_y=1,n_p
                     do p_x=1,n_p
                        p_tmp=p_tmp+1
                        pos_x_tmp(p_tmp)=pos_x(p)+(p_x-1)*delta+off
                        pos_y_tmp(p_tmp)=pos_y(p)+(p_y-1)*delta+off
                        pos_z_tmp(p_tmp)=pos_z(p)+(p_z-1)*delta+off
                        particle_mass_tmp(p_tmp)=particle_mass(p)*small
c     write(95,10)pp,p,p_tmp,
c     >                       pos_x(p),pos_y(p),pos_z(p),
c     >                       pos_x_tmp(p_tmp),pos_y_tmp(p_tmp),
c     >                       pos_z_tmp(p_tmp),
c     >                       particle_mass(p),particle_mass_tmp(p_tmp)
                     end do
                  end do
               end do
               p_tmp=p_tmp+1
               pos_x_tmp(p_tmp)=pos_x(p)
               pos_y_tmp(p_tmp)=pos_y(p)
               pos_z_tmp(p_tmp)=pos_z(p)
               particle_mass_tmp(p_tmp)=-particle_mass(p)
c               write(95,10)pp,p,p_tmp,pos_x(p),pos_y(p),pos_z(p),
c     >              pos_x_tmp(p_tmp),pos_y_tmp(p_tmp),pos_z_tmp(p_tmp),
c     >              particle_mass(p),particle_mass_tmp(p_tmp)
            end if
         end if
      end do
c     
      number_particles_tmp=p_tmp
c     
 10   format(3i8,8(1pe13.5))
      return
      end
c
