      subroutine heavies(heavy_version,grid_length,pos_x,pos_y,pos_z,
     >  pos_x_tmp,pos_y_tmp,pos_z_tmp,particle_mass,particle_mass_tmp,
     >  level,number_particles,number_particles_tmp,
     >  highest_level_group,level_min,level_max,force_max,
     >  particles_heavy_tmp,groups_maxx,particles_real)
c
      implicit none
c
      integer particles_real,groups_maxx,grid_length
      integer lev,level_max,p,number_particles,h_p,level_min
      integer highest_level_group(particles_real),p_tmp,n_p,p_x,p_y,p_z
      integer level(groups_maxx),number_particles_tmp
      integer particles_heavy_tmp
c
      logical heavy_version
c
      real log4,tmp_mass,very_small,spacing(0:8),force_max
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
      heavy_version=.false.
      log4=alog(4.0)
      tmp_mass=force_max/float(grid_length**2)
      very_small=1.0e-30
      number_particles_tmp=0
c
      print*,'heavies ',log4,tmp_mass,number_particles
      do p=1,number_particles
c        if(mod(p,100) .eq. 1) print*,' p ',p,highest_level_group(p)
        lev=level(highest_level_group(p))
        h_p=alog(tmp_mass/(particle_mass(p)+very_small))/log4
        if(h_p .lt. level_min) stop 'particle too heavy'
        heavy_version=heavy_version .or. h_p .lt. lev
      end do
c
      if(.not. heavy_version) return
c
      p_tmp=0
      do p=1,number_particles
        lev=level(highest_level_group(p))
        h_p=alog(tmp_mass/(particle_mass(p)+very_small))/log4
        if(h_p .lt. lev) then
          n_p=2**(lev-h_p)
          if(p_tmp+n_p**3+1 .gt. particles_heavy_tmp) 
     >      stop 'too many heavies'
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
              end do
            end do
          end do
          p_tmp=p_tmp+1
          pos_x_tmp(p_tmp)=pos_x(p)
          pos_y_tmp(p_tmp)=pos_y(p)
          pos_z_tmp(p_tmp)=pos_z(p)
          particle_mass_tmp(p_tmp)=-particle_mass(p)
        end if
      end do
c
      number_particles_tmp=p_tmp
c
 95   format(i8,3(1pe14.5),i8,4(1pe13.5))
      return
      end
c
