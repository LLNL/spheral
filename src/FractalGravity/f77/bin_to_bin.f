      program bin_to_bin
c     
      implicit none
c     
      integer parts_max,groups_maxx
      parameter (parts_max=550000,groups_maxx=20000)
      character*80 file1
c
      integer record_1,record_2,record_d,n,nr,number_particles,group
      integer highest_level_group(groups_maxx),mother_group(groups_maxx)
      integer groups_max
c
      real time,step_length,arad,omega_0,omega_lambda,pexp
      real pos_x(parts_max),pos_y(parts_max),pos_z(parts_max)
      real vel_x(parts_max),vel_y(parts_max),vel_z(parts_max)
      real particle_mass(parts_max)
c     
      print*,'input file'
      read(*,*)file1
      open(unit=1,file=file1,form='unformatted')
      print*,'output file'
      read(*,*)file1
      open(unit=2,file=file1,form='unformatted')
c     
      print*,'record_1,record_2,record_d'     
      read(*,*)record_1,record_2,record_d
c     
      do nr=1,record_2
        if(nr .lt. record_1) then
          read(1,end=4)time,step_length,number_particles,
     >      arad,omega_0,omega_lambda,pexp
          print*,time,step_length,number_particles,
     >      arad,omega_0,omega_lambda,pexp
          read(1)
          read(1)
          read(1)
          read(1)
          read(1)
          read(1)
          read(1)
          read(1)
          read(1)
        else
          read(1,end=4)time,step_length,number_particles,
     >      arad,omega_0,omega_lambda,pexp
          print*,time,step_length,number_particles,
     >      arad,omega_0,omega_lambda,pexp
          read(1)(pos_x(n),n=1,number_particles)
          read(1)(pos_y(n),n=1,number_particles)
          read(1)(pos_z(n),n=1,number_particles)
          read(1)(vel_x(n),n=1,number_particles)
          read(1)(vel_y(n),n=1,number_particles)
          read(1)(vel_z(n),n=1,number_particles)
          read(1)(particle_mass(n),n=1,number_particles)
          read(1)(highest_level_group(n),n=1,number_particles)
          read(1)groups_max,(mother_group(group),group=1,groups_max)
        end if
c     
        if(nr .ge. record_1 .and. mod(nr-record_1,record_d) .eq. 0)
     >    then
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
          write(2)groups_max,(mother_group(group),group=1,groups_max)
c     
        end if
      end do
c     
 4    stop
      end
