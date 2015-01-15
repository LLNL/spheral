      subroutine periodic_solver(length,number_particles,
     >     density,pot,lev,debug,
     >     groups_maxx,points_maxx,particles_maxx,tweaks,
     >     parameters_integer,parameters_real,parameters_logical)
c     
      implicit none
c     
      include 'maxx.inc'
c     
      integer groups_maxx,points_maxx,particles_maxx,tweaks
      real pot(points_maxx)
      real g_hat,green_1(grid_length_max+1),pi,grav_constant
      real density(points_maxx)
      real parameters_real(1024)
c     
      integer long(3),parameters_integer(1024)
      integer number_particles,length,point,n_x,n_y,n_z,n_x_1
      integer i,j,k,g,n,nk,lev
      logical parameters_logical(1024)
      logical debug
c     
      n(i,j,k,g)=i+((j-1)+(k-1)*g)*g
      nk(i,j,k,g)=i+((j-1)+(k-1)*g)*(g+2)
c     
      long(1)=length
      long(2)=length
      long(3)=length
c     
      pi=4.0*atan(1.0)
      grav_constant=8.0*pi/(float(length))**5
      do k=1,length+1
         green_1(k)=1.0e-30+
     >        (2.0*(sin(pi*float(k-1)/float(length))))**2
      end do
c     
      if(mod(tweaks/16,2) .eq. 1) then
         call make_rho_hat(pot,
     >        lev,length,points_maxx,
     >        parameters_integer,parameters_real,parameters_logical)
      else
         do n_z=1,length
            do n_y=1,length
               do n_x=1,length
                  point=n(n_x,n_y,n_z,length)
                  pot(point)=density(point)
               end do
            end do
         end do
c     
         call four2(pot,long,3,-1,0)
c     
      end if
c     
      if(mod(tweaks/32,2) .eq. 1) call power_spectrum(length,
     >     pot,lev,points_maxx)
c     
      pot(1)=0.0
c     
      do n_z=1,length
         do n_y=1,length
            n_x_1=0
            do n_x=1,length+2,2
               n_x_1=n_x_1+1
               g_hat=-grav_constant/
     >              (green_1(n_z)+green_1(n_y)+green_1(n_x_1))
               point=nk(n_x,n_y,n_z,length)
               pot(point)=pot(point)*g_hat
               pot(point+1)=pot(point+1)*g_hat
            end do
         end do
      end do
c     
      call four2(pot,long,3,1,-1)
c     
      return
      end

