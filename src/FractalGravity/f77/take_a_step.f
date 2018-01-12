      subroutine take_a_step(time,arad,number_particles,step_length,
     >  periodic,omega_0,omega_lambda,pexp,
     >  pos_x,pos_y,pos_z,vel_x,vel_y,vel_z,
     >  force_x,force_y,force_z,particles_max)
c     
      implicit none
c     
      integer number_particles,n,particles_max
      real pos_x(particles_max),pos_y(particles_max)
      real pos_z(particles_max),vel_x(particles_max)
      real vel_y(particles_max),vel_z(particles_max)
      real force_x(particles_max),force_y(particles_max)
      real force_z(particles_max),step_length,pexp,parad
      real omega_0,omega_lambda,time,arad,arad_half
      real v_const,f_const,h,a,om0,oml,dadt,dt,dpda,parad_half
c     
      logical periodic
c     
      h(a,om0,oml)=sqrt(om0/a**3+(1.0-om0-oml)/a**2+oml)
c
      parad=1.0
      print*,'step ',time,number_particles
      if(periodic) then
         parad=arad**pexp
         dpda=pexp*parad/arad
         parad_half=parad+step_length*0.5
         arad_half=parad_half**(1.0/pexp)
         dadt=arad_half*h(arad_half,omega_0,omega_lambda)
         v_const=1.0-2.0*step_length/arad/dpda
         f_const=step_length/arad**4/
     >        h(arad,omega_0,omega_lambda)/dpda
         dt=step_length/dadt/dpda
         do n=1,number_particles
            vel_x(n)=vel_x(n)*v_const+f_const*force_x(n)
            vel_y(n)=vel_y(n)*v_const+f_const*force_y(n)
            vel_z(n)=vel_z(n)*v_const+f_const*force_z(n)
            pos_x(n)=pos_x(n)+vel_x(n)*dt
            pos_y(n)=pos_y(n)+vel_y(n)*dt
            pos_z(n)=pos_z(n)+vel_z(n)*dt
         end do
         time=time+step_length/dadt/dpda
         parad=parad+step_length
         arad=parad**(1.0/pexp)
c     
      else
c     
         do n=1,number_particles
            if(min(force_x(n),force_y(n),force_z(n)) .gt. -1.0e29) then
               vel_x(n)=vel_x(n)+step_length*force_x(n)
               vel_y(n)=vel_y(n)+step_length*force_y(n)
               vel_z(n)=vel_z(n)+step_length*force_z(n)
               pos_x(n)=pos_x(n)+step_length*vel_x(n)
               pos_y(n)=pos_y(n)+step_length*vel_y(n)
               pos_z(n)=pos_z(n)+step_length*vel_z(n)
            end if
         end do
         time=time+step_length
      end if
c     
      return
      end
