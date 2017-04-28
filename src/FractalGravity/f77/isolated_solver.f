      subroutine isolated_solver(grid_length,density,
     >  potential_point,debug,memory_value,green,rho,
     >  groups_maxx,points_maxx,particles_maxx,tweaks)
c     
      implicit none
c     
      include 'maxx.inc'
      include 'maxx_iso.inc'
c     
      integer groups_maxx,points_maxx,particles_maxx,grid_length
      real potential_point(points_maxx)
c     
      real rho(iso_maxx)
      real green(green_length)
      real density(points_maxx)
      real y_2,z_2,r_2,grav_constant,gre
c     
      integer n_1,n_2,i,j,k,g,n,nk,tweaks,memory_value
      integer n_g_x,n_g_y,n_g_z,n_g_p,n_g_r
      integer grid_length_2,grid_length_4,n_inv,grid_length_6
      integer n_g,n_x_inv,n_y_inv,n_z_inv,grid_length_0
      integer point,n_x,n_y,n_z,grid_length_1
      integer point_0,point_x,point_y,point_z,point_g
      integer n_y_1,n_z_1,long(3)
c
      logical first_time,debug,wrap_y,wrap_z
      n(i,j,k,g)=i+((j-1)+(k-1)*g)*g
      nk(i,j,k,g)=i+((j-1)+(k-1)*g)*(g+2)
c     
      data first_time/.true./
      data grid_length_0/0/
c     
      if((2*grid_length)**2*(2*grid_length+2) .gt. iso_maxx) then
        memory_value=1
        print*,'case b memory value= ',memory_value,iso_maxx,
     >    (2*grid_length)**2*(2*grid_length+2)
        return
      end if
c
      long(1)=grid_length*2
      long(2)=grid_length*2
      long(3)=grid_length*2
c
      grav_constant=1.0/4.0/(float(grid_length))**5
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
c     
        call set_constant_real(rho,1,grid_length_6,1,0.0)
c     
        do n_z=1,grid_length_1
          z_2=(float(n_z-1))**2+0.25
          do n_y=1,grid_length_1
            y_2=(float(n_y-1))**2
            do n_x=1,grid_length_1
              r_2=z_2+y_2+(float(n_x-1))**2
              n_g=n(n_x,n_y,n_z,grid_length_2)
              rho(n_g)= -grav_constant/sqrt(r_2)
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
              rho(n_2)=rho(n_1)
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
              rho(n_2)=rho(n_1)
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
              rho(n_2)=rho(n_1)
            end do
          end do
        end do
c
        call four2(rho,long,3,-1,0)
c     
        point_g=0
        do n_z=1,grid_length_1
          do n_y=1,grid_length_1
            do n_x=1,grid_length_2+2,2
              point_g=point_g+1
              point=nk(n_x,n_y,n_z,grid_length_2)
              green(point_g)=rho(point)
            end do
          end do
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
      call four2(rho,long,3,-1,0)
c

      point_g=0
      do n_z=1,grid_length_1
        n_z_1=grid_length_2+2-n_z
        wrap_z=n_z .gt. 1 .and. n_z .lt. grid_length_1
        do n_y=1,grid_length_1
          n_y_1=grid_length_2+2-n_y
            wrap_y=n_y .gt. 1 .and. n_y .lt. grid_length_1
          do n_x=1,grid_length_2+2,2
            point_g=point_g+1
            gre=green(point_g)
c
            point=nk(n_x,n_y,n_z,grid_length_2)
            rho(point)=rho(point)*gre
            rho(point+1)=rho(point+1)*gre
c
            if(wrap_z) then
              point=nk(n_x,n_y,n_z_1,grid_length_2)
              rho(point)=rho(point)*gre
              rho(point+1)=rho(point+1)*gre
              if(wrap_y) then
                point=nk(n_x,n_y_1,n_z_1,grid_length_2)
                rho(point)=rho(point)*gre
                rho(point+1)=rho(point+1)*gre
              end if
            end if
c
            if(wrap_y) then
              point=nk(n_x,n_y_1,n_z,grid_length_2)
              rho(point)=rho(point)*gre
              rho(point+1)=rho(point+1)*gre
            end if
          end do
        end do
      end do
c     
      call four2(rho,long,3,1,-1)
c     
      do n_g_z=1,grid_length_1
        do n_g_y=1,grid_length_1
          do n_g_x=1,grid_length_1
            n_g_p=n(n_g_x,n_g_y,n_g_z,grid_length_1)
            n_g_r=n(n_g_x,n_g_y,n_g_z,grid_length_2)
            potential_point(n_g_p)=rho(n_g_r)
          end do
        end do
      end do
c     
      return
      end

