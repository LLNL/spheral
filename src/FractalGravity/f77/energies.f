      subroutine energies(time,arad,n_step,number_particles,potential,
     >  periodic,omega_0,omega_lambda,pexp,
     >  pos_x,pos_y,pos_z,vel_x,vel_y,vel_z,
     >  particle_mass,force_x,force_y,force_z,
     >  step_length,particles_max)
c     
      implicit none
c     
      integer number_particles,n_step,particles_max,p
      integer sum_0
      real potential(particles_max),vel_x(particles_max)
      real vel_y(particles_max),vel_z(particles_max)
      real particle_mass(particles_max)
      real total_energy,potential_energy,kinetic_energy
      real kinetic_energy_old
      real sum_x,sum_y,sum_z,d_x,d_y,d_z
      real pos_x(particles_max),pos_y(particles_max)
      real pos_z(particles_max)
      real force_x(particles_max),force_y(particles_max)
      real force_z(particles_max),step_length
      real omega_0,omega_lambda,time,arad,udda,a
      real sum_e_x,sum_e_y,sum_e_z,arad_up,arad_down
      real potential_energy_old
      real sum_p_x,sum_p_y,sum_p_z,sum_m,sum_pot,da
      real ang_mom_x,ang_mom_y,ang_mom_z,e_k,e_p
      real b,c,d,pexp,dp,parad_up,parad_down,parad
c     
      logical periodic,check_ok
c     
      check_ok(a,b,c,d)=min(a,b,c) .gt. d
c
      save udda,kinetic_energy_old,potential_energy_old
c
      sum_m=0.0
      sum_x=0
      sum_y=0
      sum_z=0
      sum_e_x=0
      sum_e_y=0
      sum_e_z=0
      sum_p_x=0
      sum_p_y=0
      sum_p_z=0
      sum_pot=0.0
      ang_mom_x=0.0
      ang_mom_y=0.0
      ang_mom_z=0.0
c
      sum_0=0
c
      do p=1,number_particles
         if(check_ok(force_x(p),force_y(p),force_z(p),-1.0e10)) then
            sum_0=sum_0+1
            sum_m=sum_m+particle_mass(p)
            sum_x=sum_x+pos_x(p)*particle_mass(p)
            sum_y=sum_y+pos_y(p)*particle_mass(p)
            sum_z=sum_z+pos_z(p)*particle_mass(p)
            sum_p_x=sum_p_x+vel_x(p)*particle_mass(p)
            sum_p_y=sum_p_y+vel_y(p)*particle_mass(p)
            sum_p_z=sum_p_z+vel_z(p)*particle_mass(p)
            sum_e_x=sum_e_x+vel_x(p)**2*particle_mass(p)
            sum_e_y=sum_e_y+vel_y(p)**2*particle_mass(p)
            sum_e_z=sum_e_z+vel_z(p)**2*particle_mass(p)
            sum_pot=sum_pot+potential(p)*particle_mass(p)
         end if
      end do
c
      print*,'particles inside ',n_step,sum_0
      sum_x=sum_x/sum_m
      sum_y=sum_y/sum_m
      sum_z=sum_z/sum_m
c
      do p=1,number_particles
         if(check_ok(force_x(p),force_y(p),force_z(p),-1.0e10) )then
            d_x=pos_x(p)-sum_x
            d_y=pos_y(p)-sum_y
            d_z=pos_z(p)-sum_z
            ang_mom_x=ang_mom_x+(d_y*vel_z(p)-d_z*vel_y(p))*
     >           particle_mass(p)
            ang_mom_y=ang_mom_y+(d_z*vel_x(p)-d_x*vel_z(p))*
     >           particle_mass(p)
            ang_mom_z=ang_mom_z+(d_x*vel_y(p)-d_y*vel_x(p))*
     >           particle_mass(p)
         end if
      end do
c
      potential_energy=0.5*sum_pot
      kinetic_energy=0.5*(sum_e_x+sum_e_y+sum_e_z)
c
      parad=arad**pexp
c
      if(n_step .eq. 1 .and. periodic) then
         dp=step_length
         parad_down=parad-0.5*dp
         arad_down=parad_down**(1.0/pexp)
         potential_energy_old=potential_energy*arad_down**2
         udda=-potential_energy_old*arad_down/3.0
      end if
c     
      if(periodic) then
         dp=step_length
         parad_down=parad-0.5*dp
         parad_up=parad+0.5*dp
         arad_down=parad_down**(1.0/pexp)
         arad_up=parad_up**(1.0/pexp)
         da=arad-(parad-dp)**(1.0/pexp)
c
         e_p=potential_energy*arad
         e_k=(kinetic_energy_old*arad_down**4+kinetic_energy*arad_up**4)
     >        /2.0
         udda=udda-0.5*da*(potential_energy+potential_energy_old)
         total_energy=e_p+e_k+udda
         if(n_step .gt. 0)write(3,4)n_step,time,arad,
     >        -e_p,e_k,udda,total_energy,
     >        sum_p_x,sum_p_y,sum_p_z,
     >        ang_mom_x,ang_mom_y,ang_mom_z
      else
         total_energy=0.5*(kinetic_energy_old+kinetic_energy)+
     >        potential_energy
         if(n_step .gt. 0) write(3,4)n_step,time,
     >        potential_energy,(kinetic_energy+kinetic_energy_old)*0.5,
     >        total_energy,sum_p_x,sum_p_y,sum_p_z,
     >        ang_mom_x,ang_mom_y,ang_mom_z
 4       format(i5,12(1pe11.3))
      end if
      kinetic_energy_old=kinetic_energy
      potential_energy_old=potential_energy
      return
      end
