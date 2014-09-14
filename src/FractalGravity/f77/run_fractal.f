      program run_fractal
c     
      implicit none
c     
      include 'maxes.inc'
c
      character*80 input_data_file,output_data_file
      character*80 output_energy_file
c     
      integer number_steps,output_step,minimum_number,i,j,padding
      integer level_max,moat_0,random_offset,memory_value,level_min
      integer grid_length,n_step,number_particles,n,tweaks,tweaks_0
      integer parameters_integer(1024),number_masks,level_mask(5000)
      integer maxits_sor,mother_group(groups_max),g_max
      integer highest_level_group(particles_max_0),group,groups_maxx
c     
      logical parameters_logical(1024)
      logical periodic,debug
c     
      real parameters_real(1024),pos_mask(20000)
      real pos_x(particles_max_0),pos_y(particles_max_0)
      real pos_z(particles_max_0),potential(particles_max_0)
      real force_x(particles_max_0),force_y(particles_max_0)
      real force_z(particles_max_0)
      real vel_x(particles_max_0),vel_y(particles_max_0)
      real vel_z(particles_max_0),particle_mass(particles_max_0)
      real time,step_length,arad,omega_0,omega_lambda,force_max
      real epsilon_sor,pexp
c
      data parameters_integer/1024*0/
      data parameters_real/1024*0.0/
      data parameters_logical/1024*.false./
c
      print*,'input data file'
      read(*,*)input_data_file
      print*,'output data file'
      read(*,*)output_data_file
      print*,'output energy file'
      read(*,*)output_energy_file
      print*,'number_steps,output_step'
      read(*,*)number_steps,output_step
      print*,'periodic'
      read(*,*)periodic
      print*,'minimum_number'
      read(*,*)minimum_number
      print*,'grid_length'
      read(*,*)grid_length
      print*,'level_min,level_max'
      read(*,*)level_min,level_max
      print*,'padding'
      read(*,*)padding
      print*,'moat_0'
      read(*,*)moat_0
      print*,'maximum force per particle'
      read(*,*)force_max
      print*,'random_offset'
      read(*,*)random_offset
      print*,'debug'
      read(*,*)debug
      print*,'tweaks'
      read(*,*)tweaks_0
      print*,'maxits_sor,epsilon_sor'
      read(*,*)maxits_sor,epsilon_sor
      print*,'number_masks'
      read(*,*)number_masks
      if(number_masks .gt. 0) then
         do i=1,number_masks
            read(*,*)level_mask(2*i-1),level_mask(2*i),
     >           (pos_mask(j),j=(i-1)*6+1,i*6)
         end do
      end if
c
      parameters_integer(1)=maxits_sor
      parameters_integer(4)=level_min
      parameters_real(5)=force_max
      parameters_real(8)=epsilon_sor
c
      open(unit=1,file=input_data_file,form='unformatted')
      open(unit=2,file=output_data_file,form='unformatted')
      open(unit=3,file=output_energy_file,form='formatted')
      open(unit=39,file='daughter_time.dat',form='formatted')
      if(debug) open(unit=41,file='debug.dat',form='formatted')
      open(unit=54,file='timing.dat',form='formatted')
      open(unit=68,file='level_data.dat',form='formatted')
      open(unit=93,file='sor.dat',form='formatted')
      open(unit=94,file='memory_load.dat',form='formatted')
      if(mod(tweaks/4,2) .eq. 1) then
         open(unit=97,file='groups.dat',form='formatted')
         open(unit=98,file='points.dat',form='formatted')
         open(unit=99,file='particles.dat',form='formatted')
      end if
c     
      print*,'files opened'
      do i=1,1000
        read(1,end=1)time,step_length,number_particles,
     >    arad,omega_0,omega_lambda,pexp
        read(1)(pos_x(n),n=1,number_particles)
        read(1)(pos_y(n),n=1,number_particles)
        read(1)(pos_z(n),n=1,number_particles)
        read(1)(vel_x(n),n=1,number_particles)
        read(1)(vel_y(n),n=1,number_particles)
        read(1)(vel_z(n),n=1,number_particles)
        read(1)(particle_mass(n),n=1,number_particles)
        read(1)(highest_level_group(n),n=1,number_particles)
        read(1)groups_maxx,(mother_group(group),group=1,groups_maxx)
      end do
 1    print*,'data read ',time,step_length,number_particles
c     
          print*,'printing ',n_step,time
          write(2)time,step_length,number_particles,
     >      arad,omega_0,omega_lambda,pexp
          write(2)(pos_x(n),n=1,number_particles)
          write(2)(pos_y(n),n=1,number_particles)
          write(2)(pos_z(n),n=1,number_particles)
          write(2)(vel_x(n),n=1,number_particles)
          write(2)(vel_y(n),n=1,number_particles)
          write(2)(vel_z(n),n=1,number_particles)
          write(2)(particle_mass(n),n=1,number_particles)
          write(2)(highest_level_group(n),n=1,number_particles)
          write(2)groups_maxx,(mother_group(group),group=1,groups_maxx)
c
      do n=1,number_particles
         highest_level_group(n)=-n
      end do
      do group=1,groups_max
         mother_group(group)=-group
      end do
c
      n_step=0
      print*,'energies'
      call energies(time,arad,n_step,number_particles,potential,
     >  periodic,omega_0,omega_lambda,pexp,
     >  pos_x,pos_y,pos_z,vel_x,vel_y,vel_z,
     >  particle_mass,force_x,force_y,force_z,
     >  step_length,particles_max)
c     
      print*,'data read ',time,step_length,number_particles
c     
      do n_step=1,number_steps
        print*,'fractal_gravity ',n_step,time,step_length,arad
        tweaks=tweaks_0
c
        call fractal_gravity(number_particles,grid_length,periodic,
     >    minimum_number,level_max,pos_x,pos_y,pos_z,particle_mass,
     >    potential,force_x,force_y,force_z,
     >    highest_level_group,mother_group,
     >    padding,number_masks,level_mask,pos_mask,
     >    moat_0,random_offset,debug,tweaks,memory_value,
     >    parameters_integer,parameters_real,parameters_logical)
c     
        print*,'take a step ',n_step,time,step_length,arad
        if(memory_value .eq. 0) then
          call take_a_step(time,arad,number_particles,step_length,
     >      periodic,omega_0,omega_lambda,pexp,
     >      pos_x,pos_y,pos_z,vel_x,vel_y,vel_z,
     >      force_x,force_y,force_z,particles_max)
c     
          print*,'energies ',n_step,time,arad
          call energies(time,arad,n_step,number_particles,potential,
     >      periodic,omega_0,omega_lambda,pexp,
     >      pos_x,pos_y,pos_z,vel_x,vel_y,vel_z,
     >      particle_mass,force_x,force_y,force_z,
     >      step_length,particles_max)
        end if
c     
        if((mod(n_step-1,output_step) .eq. 0 .and. n_step .gt. 1) .or. 
     >    n_step .eq. number_steps .or. memory_value .ne. 0) then
c
           g_max=1
           do group=1,groups_max
              if(mother_group(group) .gt. 0) g_max=g_max+1
           end do
c
          print*,'printing ',n_step,time,g_max
          write(2)time,step_length,number_particles,
     >      arad,omega_0,omega_lambda,pexp
          write(2)(pos_x(n),n=1,number_particles)
          write(2)(pos_y(n),n=1,number_particles)
          write(2)(pos_z(n),n=1,number_particles)
          write(2)(vel_x(n),n=1,number_particles)
          write(2)(vel_y(n),n=1,number_particles)
          write(2)(vel_z(n),n=1,number_particles)
          write(2)(particle_mass(n),n=1,number_particles)
          write(2)(highest_level_group(n),n=1,number_particles)
          write(2)g_max,(mother_group(group),group=1,g_max)
c     
          if(memory_value .ne. 0) then
            print*,'memory_value= ',memory_value
            stop 'memory value'
          end if
        end if
      end do
c     
      stop 'it is all over'
      end
c
