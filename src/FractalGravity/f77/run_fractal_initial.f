      program run_fractal_initial
c     
      implicit none
c     
      integer particles_max,groups_max
      parameter (particles_max=800000)
      parameter (groups_max=10000)
c     
      logical debug,periodic,gauss_cut,ok
      logical parameters_logical(1024)
      logical dist_p
      character*80 output_data_file,input_level_file
c     
      real parameters_real(1024),ratio_rms_rho,ratio_rms_force
      real omega_b,h,box_length,scaling,a1,a2,alpha,epsilon_sor
      real omega_0,omega_lambda,offset_x,offset_y,offset_z,varv
      real delta,mass,rms_rho,rms_vel,rms_force_obs,rms_force,wave_0
      real pos_x(particles_max),pos_y(particles_max)
      real pos_z(particles_max),potential(particles_max)
      real pos_x_tmp(particles_max),pos_y_tmp(particles_max)
      real pos_z_tmp(particles_max)
      real vel_x(particles_max),vel_y(particles_max)
      real vel_z(particles_max),particle_mass(particles_max)
      real force_x(particles_max),force_y(particles_max)
      real force_z(particles_max),arad,time,step_length,pi,rms
      real force_x_tmp(particles_max),force_y_tmp(particles_max)
      real force_z_tmp(particles_max),delta_l(0:10),left_l(0:10)
      real var_force,ratio_rho,ratio_force,pos_mask(20000),hubble,omega
      real extra_factor,da,om,force_max,growth,amplitude,redshift,lambda
      real omega_z,lambda_z,h_z,pexp
      real scale,power_slope,cut_off,aa,bb,cc
      real conv_f_p,conv_f_v,rms_rho_obs,rr,ratio
      real conv,var_rho,age_of_the_universe,dp
      real rms_density,rms_velocity,rms_vel_obs,dist_1
c
      integer p_x,p_y,p_z,normalization,idum(0:8),spectrum_number,mm,m
      integer p_up_x,p_down_x,p_y_up,p_y_down,p_z_up,p_z_down
      integer parameters_integer(1024),maxits_sor
      integer grid_length,particles,number_particles,memory_value,tweaks
      integer n_x,n_y,n_z,p,padding,moat_0,random_offset
      integer minimum_number,level_max,n,nn,j,k,l,particle,i
      integer mother_group(groups_max),lev,level_max_zel
      integer new_l(0:10),np,level_min,groups_maxx
      integer highest_level_group(particles_max),level(particles_max)
      integer number_masks,level_mask(4000),max_lev_a,max_lev_b,group
c     
      equivalence (vel_x(1),pos_x_tmp(1))
      equivalence (vel_y(1),pos_y_tmp(1))
      equivalence (vel_z(1),pos_z_tmp(1))
      equivalence (vel_x(1),force_x_tmp(1))
      equivalence (vel_y(1),force_y_tmp(1))
      equivalence (vel_z(1),force_z_tmp(1))
c     
      dist_1(aa,bb)=min(abs(aa-bb),abs(aa-bb+1.0),abs(aa-bb-1.0))
      dist_p(aa,bb,cc)=dist_1(aa,bb) .le. cc
c
      nn(i,j,k,l)=i+((j-1)+(k-1)*l)*l
c
      hubble(arad,omega_0,omega_lambda)=sqrt(omega_0/arad**3+
     >     (1.0-omega_0-omega_lambda)/arad**2+omega_lambda)
c
      omega(arad,omega_0,omega_lambda)=
     >     omega_0/(arad**3*hubble(arad,omega_0,omega_lambda)**2)
c
      lambda(arad,omega_0,omega_lambda)=
     >     omega_lambda/hubble(arad,omega_0,omega_lambda)**2
c
      growth(arad,omega_0,omega_lambda)=
     >     2.5*omega(arad,omega_0,omega_lambda)*arad/
     >     (omega(arad,omega_0,omega_lambda)**(4.0/7.0)-
     >     lambda(arad,omega_0,omega_lambda)+
     >     (1.0+0.5*omega(arad,omega_0,omega_lambda))*
     >     (1.0+lambda(arad,omega_0,omega_lambda/70.0)))
c
      amplitude(arad,omega_0,omega_lambda)=
     >     growth(arad,omega_0,omega_lambda)/
     >     growth(1.0,omega_0,omega_lambda)
c     
      data padding,moat_0,periodic,minimum_number/0,0,.true.,8/
      data random_offset/0/
c     
      data parameters_integer/1024*0/
      data parameters_real/1024*0.0/
      data parameters_logical/1024*.false./
c
      print*,'output file'
      read(*,*)output_data_file
c     
      print*,'input level file'
      read(*,*)input_level_file
c     
      print*,'delta t, pexp'
      read(*,*)step_length,pexp
      print*,step_length,pexp
c     
      print*,'omega_0,omega_lambda,redshift'
      read(*,*)omega_0,omega_lambda,redshift
      print*,omega_0,omega_lambda,redshift
c     
      print*,'omega_b,h,box length'
      read(*,*)omega_b,h,box_length
      print*,omega_b,h,box_length
c     
      print*,'spectrum number'
      read(*,*)spectrum_number
      print*,spectrum_number
c     
      print*,'power spectrum slope'
      read(*,*)power_slope
      print*,power_slope
c     
      print*,'cut_off'
      read(*,*)cut_off
      print*,cut_off
c     
      print*,'normalization,scale,rms'
      read(*,*)normalization,scale,rms
      print*,normalization,scale,rms
c     
      print*,'grid length'
      read(*,*)grid_length
      print*,grid_length
c     
      print*,'level_min,level_max,level_max_zel'
      read(*,*)level_min,level_max,level_max_zel
      print*,level_min,level_max,level_max_zel
c     
      level_max_zel=min(level_max_zel,level_max)
c
      print*,'particles across'
      read(*,*)particles
      print*,particles
c     
      print*,'maxits_sor,epsilon_sor'
      read(*,*)maxits_sor,epsilon_sor
      print*,maxits_sor,epsilon_sor
c
      print*,'maximum force per particle'
      read(*,*)force_max
      print*,force_max
c
      print*,'offset in grid units'
      read(*,*)offset_x,offset_y,offset_z
      print*,offset_x,offset_y,offset_z
c     
      print*,'debug'
      read(*,*)debug
      print*,debug
c     
      print*,'random integers for power spectrum'
      read(*,*)(idum(lev),lev=0,max(0,level_max_zel))
      print*,(idum(lev),lev=0,max(0,level_max_zel))
c     
      print*,'gaussian cut'
      read(*,*)gauss_cut
      print*,gauss_cut
c     
      print*,'number_masks'
      read(*,*)number_masks
      if(number_masks .gt. 0) then
        do i=1,number_masks
          read(*,*)level_mask(2*i-1),level_mask(2*i),
     >          (pos_mask(j),j=(i-1)*6+1,i*6)
        end do
      else
           number_masks=1
           level_mask(1)=0
           level_mask(2)=0
           pos_mask(1)=0.5
           pos_mask(2)=0.5
           pos_mask(3)=0.5
           pos_mask(4)=2.0
           pos_mask(5)=2.0
           pos_mask(6)=2.0
      end if
c     
      print*,'here 1 '
      open(unit=2,file=output_data_file,form='unformatted')
      open(unit=54,file='timing.dat',form='formatted')
c     
      number_particles=particles**3
      do particle=1,number_particles
         level(particle)=0
      end do
c
      print*,'here 2 '
      groups_maxx=groups_max
c
      h_z=hubble(1.0/(1.0+redshift),omega_0,omega_lambda)
      omega_z=omega(1.0/(1.0+redshift),omega_0,omega_lambda)
      lambda_z=lambda(1.0/(1.0+redshift),omega_0,omega_lambda)
      rms=rms*amplitude(1.0/(1.0+redshift),omega_0,omega_lambda)
c
      print*,'omega ',omega_0,omega_lambda,redshift,h_z,omega_z,lambda_z
c
      open(unit=1,file=input_level_file,form='unformatted',status='old',
     >     err=3)
      read(1)
      read(1)
      read(1)
      read(1)
      read(1)
      read(1)
      read(1)
      read(1)
        read(1)(highest_level_group(n),n=1,number_particles)
        read(1)groups_maxx,(mother_group(group),group=1,groups_maxx)
c
        do n=1,number_particles
           level(n)=-1
           group=highest_level_group(n)
           do while(group .gt. 0)
              level(n)=level(n)+1
              group=mother_group(group)
           end do
        end do
c
 3    rms_rho=0.0
      rms_vel=0.0
      if(normalization .eq. 0) then
         rms_rho=rms
      else
         rms_vel=rms
      end if
c     
      delta=1.0/float(particles)
c
      do lev=0,10
        delta_l(lev)=2.0**(-lev)*delta
        left_l(lev)=-delta_l(lev)*float(2**lev-1)/2.0
        new_l(lev)=2**lev
c        write(31,*)lev,delta_l(lev),left_l(lev),new_l(lev)
      end do
c     
      
      if(spectrum_number .eq. 0) then
        scaling=1.0
      else if(spectrum_number .eq. 1) then
c     
        a1=(46.9*omega_0*h*h)**0.67*
     >     (1.0+(32.1*omega_0*h*h)**(-0.532))
        a2=(12.0*omega_0*h*h)**0.424*
     >     (1.0+(45.0*omega_0*h*h)**(-0.582))
        alpha=a1**(-omega_b/omega_0)*a2**(-(omega_b/omega_0)**3)
        scaling=box_length*omega_0*h**2*sqrt(alpha)*
     >     (1.0-omega_b/omega_0)**0.6
        print*,'cdm ',omega_0,omega_lambda,omega_b,h,a1,a2,
     >     alpha,scaling
      else
        print*,'wrong spectrum'
        stop 'die you scumsucking pig'
      end if
c     
      pi=4.0*atan(1.0)
      time=age_of_the_universe(omega_z,lambda_z)
      arad=1.0
      mass=3.0/8.0/pi*omega_0/float(number_particles)
c     
      n=0
      do n_z=1,particles
        do n_y=1,particles
          do n_x=1,particles
            n=n+1
            pos_x(n)=(float(n_x-1)+offset_x)*delta
            pos_y(n)=(float(n_y-1)+offset_y)*delta
            pos_z(n)=(float(n_z-1)+offset_z)*delta
            particle_mass(n)=mass
          end do
        end do
      end do
c     
      memory_value=0
      tweaks=48
      wave_0=1.0/scale
      parameters_real(1)=power_slope
      parameters_real(2)=cut_off
      parameters_real(3)=wave_0
      parameters_real(4)=scaling
      parameters_real(5)=force_max
      parameters_real(8)=epsilon_sor
      parameters_logical(1)=gauss_cut
      parameters_integer(1)=maxits_sor
      parameters_integer(2)=spectrum_number
      parameters_integer(3)=level_max_zel
      parameters_integer(4)=level_min
c
      do lev=0,max(0,level_max_zel)
         parameters_integer(100+lev)=idum(lev)
      end do
c     
      print*,'input to fractal ',number_particles,grid_length,periodic,
     >   minimum_number,level_max,padding,time,
     >   moat_0,random_offset,debug,tweaks,memory_value
c     
      do n=1,number_particles
         highest_level_group(n)=-n
      end do
      do group=1,groups_max
         mother_group(group)=-group
      end do
c
      print*,'going in 1'
      call fractal_gravity(number_particles,grid_length,periodic,
     >   minimum_number,0,pos_x,pos_y,pos_z,particle_mass,
     >   potential,force_x,force_y,force_z,
     >   highest_level_group,mother_group,
     >   padding,number_masks,level_mask,pos_mask,
     >   moat_0,random_offset,debug,tweaks,memory_value,
     >   parameters_integer,parameters_real,parameters_logical)
c     
      print*,'coming out 1'
      if(memory_value .ne. 0) stop 'not enough memory 1'
c
      ratio_rms_rho=parameters_real(7)
      ratio_rms_force=parameters_real(6)
c     
      conv=float(grid_length)/2.0/4.0/pi
c     
      var_rho=0.0
      var_force=0.0
      do p_z=1,particles
        do p_y=1,particles
          do p_x=1,particles
            p_z_up=mod(p_z,particles)+1
            p_z_down=mod(p_z-2+particles,particles)+1
            p_z_up=nn(p_x,p_y,p_z_up,particles)
            p_z_down=nn(p_x,p_y,p_z_down,particles)
            p_y_up=mod(p_y,particles)+1
            p_y_down=mod(p_y-2+particles,particles)+1
            p_y_up=nn(p_x,p_y_up,p_z,particles)
            p_y_down=nn(p_x,p_y_down,p_z,particles)
            p_up_x=mod(p_x,particles)+1
            p_down_x=mod(p_x-2+particles,particles)+1
            p_up_x=nn(p_up_x,p_y,p_z,particles)
            p_down_x=nn(p_down_x,p_y,p_z,particles)
            p=nn(p_x,p_y,p_z,particles)
            rr=(force_x(p_up_x)-force_x(p_down_x)+
     >         force_y(p_y_up)-force_y(p_y_down)+
     >         force_z(p_z_up)-force_z(p_z_down))*conv
            var_rho=var_rho+rr**2
            var_force=var_force+force_x(p)**2+force_y(p)**2+
     >         force_z(p)**2
          end do
        end do
      end do
c     
      rms_rho_obs=sqrt(var_rho/float(number_particles))*
     >   8.0*pi/3.0/omega_0*ratio_rms_rho
      rms_force_obs=sqrt(var_force/float(3*number_particles))*
     >   ratio_rms_force
      rms_vel_obs=rms_force_obs/(1.5*omega_z**0.43)
c     
      ratio_rho=rms_rho/rms_rho_obs
      rms_force=rms_vel*1.5*omega_z**0.43
      ratio_force=rms_force/rms_force_obs
      if(normalization .eq. 0) then
        ratio=ratio_rho
      else
        ratio=ratio_force
      end if
c     
      rms_density=rms_rho_obs*ratio
      rms_velocity=rms_vel_obs*ratio
      print*,'rms rho, rms vel, scale ',rms_density,rms_velocity,scale
c     
      conv=2.0/3.0/omega_z
c     
      do p=1,number_particles
        pos_x_tmp(p)=pos_x(p)
        pos_y_tmp(p)=pos_y(p)
        pos_z_tmp(p)=pos_z(p)
      end do
c     
      ok=.false.
      n=0
      np=0
      do n_z=1,particles
        do n_y=1,particles
          do n_x=1,particles
            np=np+1
            lev=level(np)
            max_lev_a=-1
            max_lev_b=-1
c     
            if(number_masks .gt. 0) then
               mm=-5
               do m=1,number_masks
                  mm=mm+6
                  if(level_mask(2*m) .eq. 0) then
                     ok=dist_p(pos_x_tmp(np),pos_mask(mm),
     >                    pos_mask(mm+3))
                     ok=ok .and. dist_p(pos_y_tmp(np),pos_mask(mm+1),
     >                    pos_mask(mm+4))
                     ok=ok .and. dist_p(pos_z_tmp(np),pos_mask(mm+2),
     >                    pos_mask(mm+5))
                  else if(level_mask(2*m) .eq. 1) then
                     ok=(dist_1(pos_x_tmp(np),pos_mask(mm))/
     >                    pos_mask(mm+3))**2+
     >                    (dist_1(pos_y_tmp(np),pos_mask(mm+1))/
     >                    pos_mask(mm+4))**2+
     >                    (dist_1(pos_z_tmp(np),pos_mask(mm+2))/
     >                    pos_mask(mm+5))**2 .le. 1.0
                  else
                     stop 'baad maasking'
               end if
c
                  if(ok) then
                     if(level_mask(2*m-1) .ge. 0) then
                        max_lev_a=max(max_lev_a,level_mask(2*m-1))
                     else
                        max_lev_b=max(max_lev_b,-level_mask(2*m-1))
                     end if
c                     write(31,*)'joke ',np,lev,max_lev_a,max_lev_b
                  end if
               end do
            end if
c     
            if(max_lev_a .ge. 0) lev=min(lev,max_lev_a)
            if(max_lev_b .ge. 0) lev=max(lev,max_lev_b)
c            write(31,*)np,n,lev,pos_x(np),pos_y(np),pos_z(np)
            do k=1,new_l(lev)
               do j=1,new_l(lev)
                  do i=1,new_l(lev)
                     n=n+1
                     if(n .gt. particles_max) then
                        print*,'np= ',np,n,lev,new_l(lev),
     >                       pos_x_tmp(np),
     >                       pos_y_tmp(np),
     >                       pos_z_tmp(np)
                        stop 'not enough particle space'
                     end if
                     pos_x(n)=pos_x_tmp(np)+
     >                    left_l(lev)+delta_l(lev)*float(i-1)
                     pos_y(n)=pos_y_tmp(np)+
     >                    left_l(lev)+delta_l(lev)*float(j-1)
                     pos_z(n)=pos_z_tmp(np)+
     >                    left_l(lev)+delta_l(lev)*float(k-1)
                     particle_mass(n)=mass/8.0**lev
                     pos_x(n)=mod(pos_x(n)+1.0,1.0)
                     pos_y(n)=mod(pos_y(n)+1.0,1.0)
                     pos_z(n)=mod(pos_z(n)+1.0,1.0)
c                     write(31,31)np,n,i,j,k,lev,level(np)
c     >                    pos_x_tmp(np),pos_y_tmp(np),pos_z_tmp(np),
c     >                    pos_x(n),pos_y(n),pos_z(n),
c     >                    particle_mass(n),1.0/delta_l(lev)
 31                  format(2i8,5i2,8(1pe13.5))
                  end do
               end do
            end do
         end do
      end do
      end do
c     
      number_particles=n
      memory_value=0
      tweaks=48
      wave_0=1.0/scale
      parameters_real(1)=power_slope
      parameters_real(2)=cut_off
      parameters_real(3)=wave_0
      parameters_real(4)=scaling
      parameters_logical(1)=gauss_cut
      parameters_integer(2)=spectrum_number
c     
      write(92)time,step_length,number_particles,
     >      arad,omega_z,lambda_z,pexp
      write(92)(pos_x(n),n=1,number_particles)
      write(92)(pos_y(n),n=1,number_particles)
      write(92)(pos_z(n),n=1,number_particles)
      write(92)(0.0,n=1,number_particles)
      write(92)(0.0,n=1,number_particles)
      write(92)(0.0,n=1,number_particles)
      write(92)(particle_mass(n),n=1,number_particles)
      write(92)(highest_level_group(n),n=1,number_particles)
      write(92)groups_max,(mother_group(group),group=1,groups_max)
c
      do lev=0,max(0,level_max_zel)
         parameters_integer(100+lev)=idum(lev)
      end do
c     
      do n=1,number_particles
         highest_level_group(n)=-n
      end do
      do group=1,groups_max
         mother_group(group)=-group
      end do
c
      call fractal_gravity(number_particles,grid_length,periodic,
     >   minimum_number,level_max,pos_x,pos_y,pos_z,particle_mass,
     >   potential,force_x,force_y,force_z,
     >   highest_level_group,mother_group,
     >   padding,number_masks,level_mask,pos_mask,
     >   moat_0,random_offset,debug,tweaks,memory_value,
     >   parameters_integer,parameters_real,parameters_logical)
c
      if(memory_value .ne. 0) stop 'not enough memory 2'
c      
      conv_f_p=ratio*2.0/3.0/omega_z
c     
      do p=1,number_particles
        force_x_tmp(p)=force_x(p)*ratio
        force_y_tmp(p)=force_y(p)*ratio
        force_z_tmp(p)=force_z(p)*ratio
c        write(47,37)p,pos_x(p),pos_y(p),pos_z(p),
c     >       force_x_tmp(p),force_y_tmp(p),force_z_tmp(p),
c     >       particle_mass(p)
c 37     format(i8,7(1pe13.5))
        pos_x(p)=mod(pos_x(p)+conv_f_p*force_x(p)+1.0,1.0)
        pos_y(p)=mod(pos_y(p)+conv_f_p*force_y(p)+1.0,1.0)
        pos_z(p)=mod(pos_z(p)+conv_f_p*force_z(p)+1.0,1.0)
      end do
c     
      memory_value=0
      tweaks=0
c     
      do n=1,number_particles
         highest_level_group(n)=-n
      end do
      do group=1,groups_max
         mother_group(group)=-group
      end do
c
      call fractal_gravity(number_particles,grid_length,periodic,
     >   minimum_number,level_max,pos_x,pos_y,pos_z,particle_mass,
     >   potential,force_x,force_y,force_z,
     >   highest_level_group,mother_group,
     >   padding,number_masks,level_mask,pos_mask,
     >   moat_0,random_offset,debug,tweaks,memory_value,
     >   parameters_integer,parameters_real,parameters_logical)
c     
      if(memory_value .ne. 0) stop 'not enough memory 3'
c
      conv_f_p=2.0/3.0/omega_z
      conv_f_v=2.0/3.0/omega_z**0.43
c
      dp=-0.5*step_length
      da=(1.0+dp)**(1.0/pexp)-1.0
c      da=-0.5*step_length
      om=omega(1.0+da,omega_z,lambda_z)
      extra_factor=(om/omega_z)**0.57*(1.0+da)*
     >     hubble(1.0+da,omega_z,lambda_z)
c     
      varv=0.0
      do p=1,number_particles
        pos_x(p)=mod(pos_x(p)+conv_f_p*(force_x(p)-force_x_tmp(p))+1.0,
     >     1.0)
        pos_y(p)=mod(pos_y(p)+conv_f_p*(force_y(p)-force_y_tmp(p))+1.0,
     >     1.0)
        pos_z(p)=mod(pos_z(p)+conv_f_p*(force_z(p)-force_z_tmp(p))+1.0,
     >     1.0)
        vel_x(p)=conv_f_v*force_x(p)*extra_factor
        vel_y(p)=conv_f_v*force_y(p)*extra_factor
        vel_z(p)=conv_f_v*force_z(p)*extra_factor
        varv=varv+vel_x(p)**2+vel_y(p)**2+vel_z(p)**2
      end do
c
      varv=sqrt(varv/float(number_particles*3))
      print*,'vel variance ',varv
c     
      write(2)time,step_length,number_particles,
     >   arad,omega_z,lambda_z,pexp
      write(2)(pos_x(n),n=1,number_particles)
      write(2)(pos_y(n),n=1,number_particles)
      write(2)(pos_z(n),n=1,number_particles)
      write(2)(vel_x(n),n=1,number_particles)
      write(2)(vel_y(n),n=1,number_particles)
      write(2)(vel_z(n),n=1,number_particles)
      write(2)(particle_mass(n),n=1,number_particles)
      write(2)(highest_level_group(n),n=1,number_particles)
      write(2)groups_max,(mother_group(group),group=1,groups_max)
c     
      print*,'particles ',number_particles
      stop 'it is all over'
      end
c
      real function age_of_the_universe(omega_0,omega_lambda)
c
      implicit none
c
      real omega_0,omega_lambda,x,dx
      integer n
c
      age_of_the_universe=0.0
      dx=1.0e-4
      do n=1,10000
         x=dx*float(n-1)+dx*0.5
         age_of_the_universe=age_of_the_universe+
     >        dx*x/sqrt(omega_0*x+(1.0-omega_0-omega_lambda)*x*x+
     >        omega_lambda*x**4)
      end do
c
      return
      end
c
