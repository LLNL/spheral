      subroutine isolated_solver(grid_length,density,
     >  potential_point,debug,memory_value,
     >  groups_maxx,points_maxx,particles_maxx,tweaks)
c     
      implicit none
c     
      include 'maxx.inc'
      include 'maxx_iso.inc'
c     
      integer groups_maxx,points_maxx,particles_maxx,grid_length
      real potential_point(points_maxx),a_grid
      real rho_nyq(iso*8*grid_length_max**2+1)
      real green_nyq(iso*4*grid_length_max**2+1)
c     
      real rho(iso_maxx),green((iso_maxx+1)/2)
      real density(points_maxx)
      real green_0(iso_maxx),green_0_nyq(iso*4*grid_length_max**2+1)
c     
      integer n_1,n_2,i,j,k,g,n,tweaks,memory_value
      integer n_g_x,n_g_y,n_g_z,n_g_p,n_g_r
      integer grid_length_2,grid_length_4,n_inv,grid_length_6
      integer n_g,n_x_inv,n_y_inv,n_z_inv,grid_length_0
      integer point,n_x,n_y,n_z,grid_length_1
      integer point_0,point_x,point_y,point_z,point_g
      real y_2,z_2,r_2,grav_constant
      save green,green_nyq
      equivalence (rho(1),green_0(1))
      equivalence (rho_nyq(1),green_0_nyq(1))
c     
      logical first_time,debug
      n(i,j,k,g)=i+((j-1)+(k-1)*g)*g
c     
      data first_time/.true./
      data grid_length_0/0/
c
c***  call SET_CONSTANT_REAL
c***  set a real array to a constant value
c***  call RLFT3
C***  CALL REAL FFT
c     
      if((2*grid_length)**3 .gt. iso_maxx) then
        memory_value=1
        print*,'case b memory value= ',memory_value,iso_maxx,
     >    (2*grid_length)**3
        return
      end if
c
      grav_constant=1.0/float(4*grid_length**2)
      a_grid=grid_length
      grid_length_1=grid_length+1
      grid_length_2=2*grid_length
      grid_length_4=grid_length_2**2
      grid_length_6=grid_length_2**3
c     
      first_time=first_time .or. grid_length .ne. grid_length_0
c     
      if(first_time) then
       grid_length_0=grid_length
       first_time=.false.
c       epsilon_2=1.0/a_grid**2
c     
       call set_constant_real(green_0,1,grid_length_6,1,0.0)
c     
       do n_z=1,grid_length+1
        z_2=(float(n_z-1))**2+0.25
        do n_y=1,grid_length+1
         y_2=(float(n_y-1))**2
         do n_x=1,grid_length+1
          r_2=z_2+y_2+(float(n_x-1))**2
          n_g=n(n_x,n_y,n_z,grid_length_2)
          green_0(n_g)= -grav_constant/sqrt(r_2)
         end do
        end do
       end do
c     
       n_inv=grid_length_2+2
c     
       do n_z=1,grid_length_1
        do n_y=1,grid_length_1
         do n_x=2,grid_length
          n_x_inv=n_inv-n_x
          n_1=n(n_x,n_y,n_z,grid_length_2)
          n_2=n(n_x_inv,n_y,n_z,grid_length_2)
          green_0(n_2)=green_0(n_1)
         end do
        end do
       end do
c     
       do n_z=1,grid_length_1
        do n_y=2,grid_length
         do n_x=1,grid_length_2
          n_y_inv=n_inv-n_y
          n_1=n(n_x,n_y,n_z,grid_length_2)
          n_2=n(n_x,n_y_inv,n_z,grid_length_2)
          green_0(n_2)=green_0(n_1)
         end do
        end do
       end do
c     
       do n_z=2,grid_length
        do n_y=1,grid_length_2
         do n_x=1,grid_length_2
          n_z_inv=n_inv-n_z
          n_1=n(n_x,n_y,n_z,grid_length_2)
          n_2=n(n_x,n_y,n_z_inv,grid_length_2)
          green_0(n_2)=green_0(n_1)
         end do
        end do
       end do
c     
       call rlft3(green_0,green_0_nyq,
     >   grid_length_2,grid_length_2,grid_length_2,1)
c     
      point_g=0
      do n_z=1,nyq+1
        do n_y=1,nyq+1
          do n_x=1,nyq
            point_g=point_g+1
            point=n(n_x,n_y,n_z,grid_length_2)
            green(point_g)=green_0(point)
          end do
        end do
      end do
c      do point=1,grid_length_6,2
c         point_g=point_g+1
c         green(point_g)=green_0(point)
c      end do
c     
      point_g=0
      do point=1,grid_length_4*2,2
         point_g=point_g+1
         green_nyq(point_g)=green_0_nyq(point_g)
      end do
c     
      end if
c     
      call set_constant_real(rho,1,grid_length_6,1,0.0)
c     
      do point_z=1,grid_length_1
       do point_y=1,grid_length_1
        do point_x=1,grid_length_1
         point=n(point_x,point_y,point_z,grid_length_2)
         point_0=n(point_x,point_y,point_z,grid_length_1)
         rho(point)=density(point_0)
        end do
       end do
      end do
c     
      call rlft3(rho,rho_nyq,
     >  grid_length_2,grid_length_2,grid_length_2,1)
c
      point_g=0
c     
      point_g=0
      do point=1,grid_length_6,2
         point_g=point_g+1
       rho(point)=rho(point)*green(point_g)
       rho(point+1)=rho(point+1)*green(point_g)
      end do
c     
      point_g=0
      do point=1,grid_length_4*2,2
         point_g=point_g+1
         rho_nyq(point)=rho_nyq(point)*green_nyq(point_g)
         rho_nyq(point+1)=rho_nyq(point+1)*green_nyq(point_g)
      end do
c     
      call rlft3(rho,rho_nyq,
     >  grid_length_2,grid_length_2,grid_length_2,-1)
c     
      do n_g_z=1,grid_length+1
       do n_g_y=1,grid_length+1
        do n_g_x=1,grid_length+1
         n_g_p=n(n_g_x,n_g_y,n_g_z,grid_length_1)
         n_g_r=n(n_g_x,n_g_y,n_g_z,grid_length_2)
         potential_point(n_g_p)=rho(n_g_r)
        end do
       end do
      end do
c     
      return
      end

