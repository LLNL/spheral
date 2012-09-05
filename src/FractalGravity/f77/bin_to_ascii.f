      program bin_to_ascii
c     
      implicit none
c     
      integer parts_max
      parameter (parts_max=800000)
      character*80 file1
      integer record_1,record_2,record_d,n,nr,number_particles,p
      integer highest_level_group(parts_max),mother_group(20000)
      integer groups_maxx,group
      real time,step_length,arad,omega_0,omega_lambda
      real pos_x(parts_max),pos_y(parts_max),pos_z(parts_max)
      real vel_x(parts_max),vel_y(parts_max),vel_z(parts_max)
      real particle_mass(parts_max)
      real xmin,xmax,ymin,ymax,zmin,zmax,pexp
c     
      print*,'input file'
      read(*,*)file1
      open(unit=1,file=file1,form='unformatted')
      print*,'output file'
      read(*,*)file1
      open(unit=2,file=file1)
c     
      print*,'record_1,record_2,record_d'     
      read(*,*)record_1,record_2,record_d
c     
      print*,'xmin,xmax,ymin,ymax,zmin,zmax'
      read(*,*)xmin,xmax,ymin,ymax,zmin,zmax
c     
      if(record_1 .gt. 1) then
         do nr=1,record_1-1
            read(1,end=4)time,step_length,number_particles,
     >           arad,omega_0,omega_lambda,pexp
            print*,time,step_length,number_particles,
     >           arad,omega_0,omega_lambda,pexp
            read(1)
            read(1)
            read(1)
            read(1)
            read(1)
            read(1)
            read(1)
            read(1)
            read(1)
         end do
      end if
      do nr=record_1,record_2
         read(1,end=4)time,step_length,number_particles,
     >        arad,omega_0,omega_lambda,pexp
         print*,time,step_length,number_particles,
     >        arad,omega_0,omega_lambda,pexp
         if(mod(nr-record_1,record_d) .eq. 0) then
            read(1)(pos_x(n),n=1,number_particles)
            read(1)(pos_y(n),n=1,number_particles)
            read(1)(pos_z(n),n=1,number_particles)
            read(1)(vel_x(n),n=1,number_particles)
            read(1)(vel_y(n),n=1,number_particles)
            read(1)(vel_z(n),n=1,number_particles)
            read(1)(particle_mass(n),n=1,number_particles)
           read(1)(highest_level_group(n),n=1,number_particles)
           read(1)groups_maxx,(mother_group(group),group=1,groups_maxx)
         else
            read(1)
            read(1)
            read(1)
            read(1)
            read(1)
            read(1)
            read(1)
            read(1)
            read(1)
         end if
c     
         if(mod(nr-record_1,record_d) .eq. 0) then
            write(2,2)number_particles,time,step_length,
     >           arad,omega_0,omega_lambda,pexp
c     
            write(2,3)xmin,ymin,zmin,0.0,0.0,0.0,0.0
            write(2,3)xmax,ymax,zmax,0.0,0.0,0.0,0.0
            do p=1,number_particles
               if(pos_x(p) .ge. xmin .and. pos_x(p) .le. xmax) then
                  if(pos_y(p) .ge. ymin .and. pos_y(p) .le. ymax) then
                     if(pos_z(p) .ge. zmin .and. pos_z(p) .le. zmax) 
     >                    then
                        write(2,3)pos_x(p),pos_y(p),pos_z(p),
     >                       vel_x(p),vel_y(p),vel_z(p),
     >                       particle_mass(p),highest_level_group(p)
                     end if
                  end if
               end if
            end do
         end if
      end do
c     
 2    format(i7,6(1pe11.3))
 3    format(7(1pe13.5),i6)
 4    stop
      end
