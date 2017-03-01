#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
#include "nbody.hh"
void energies(Nbody& nbody,Fractal& fractal)
{
  static ofstream FileEnergy;
  if(!FileEnergy.is_open())
    FileEnergy.open("energy.d");
  static ofstream FileMom;
  if(!FileMom.is_open())
    FileMom.open("momentum.d");
  double sum_m=0.0;
  double sum_x=0.0;
  double sum_y=0.0;
  double sum_z=0.0;
  double sum_vx=0.0;
  double sum_vy=0.0;
  double sum_vz=0.0;
  double sum_ekx=0.0;
  double sum_eky=0.0;
  double sum_ekz=0.0;
  double sum_p=0.0;
  for(int n=0;n < nbody.number_particles;++n)
    {
      double m=nbody.particle_mass[n];
      sum_m+=m;
      sum_x+=m*nbody.pos_x[n];
      sum_y+=m*nbody.pos_y[n];
      sum_z+=m*nbody.pos_z[n];
      sum_vx+=m*nbody.vel_x[n];
      sum_vy+=m*nbody.vel_y[n];
      sum_vz+=m*nbody.vel_z[n];
      sum_ekx+=m*pow(nbody.vel_x[n],2);
      sum_eky+=m*pow(nbody.vel_y[n],2);
      sum_ekz+=m*pow(nbody.vel_z[n],2);
      sum_p+=m*fractal.get_potential(n);
    }
  sum_x/=sum_m;
  sum_y/=sum_m;
  sum_z/=sum_m;
  double sum_px=0.0;
  double sum_py=0.0;
  double sum_pz=0.0;
  for(int n=0;n < nbody.number_particles;++n)
    {
      double m=nbody.particle_mass[n];
      sum_px+=m*((nbody.pos_y[n]-sum_y)*nbody.vel_z[n]-(nbody.pos_z[n]-sum_z)*nbody.vel_y[n]);
      sum_py+=m*((nbody.pos_z[n]-sum_z)*nbody.vel_x[n]-(nbody.pos_x[n]-sum_x)*nbody.vel_z[n]);
      sum_pz+=m*((nbody.pos_x[n]-sum_x)*nbody.vel_y[n]-(nbody.pos_y[n]-sum_y)*nbody.vel_x[n]);
    }
  nbody.potential_energy=0.5*sum_p;
  nbody.kinetic_energy=0.5*(sum_ekx+sum_eky+sum_ekz);
  if(nbody.periodic)
    {
      double parad=pow(nbody.arad,nbody.pexp);
      if(nbody.steps == 0)
	{
	  double arad_down=pow(parad-0.5*nbody.step_length,1.0/nbody.pexp);
	  nbody.potential_energy_old=nbody.potential_energy*arad_down*arad_down;
	  nbody.udda=-nbody.potential_energy_old*arad_down/3.0;
	}
      double parad_down=parad-0.5*nbody.step_length;
      double parad_up=parad+0.5*nbody.step_length;
      double arad_down=pow(parad_down,1.0/nbody.pexp);
      double arad_up=pow(parad_up,1.0/nbody.pexp);
      double da=arad_up-arad_down;
      double e_p=nbody.potential_energy*nbody.arad;
      double e_k=(nbody.kinetic_energy_old*pow(arad_down,4)+
		  nbody.kinetic_energy*pow(arad_up,4))/2.0;
      nbody.udda-=0.5*da*(nbody.potential_energy+nbody.potential_energy_old);
      double total_energy=e_p+e_k+nbody.udda;
      if(nbody.steps >= 0)
	{
	  FileEnergy << nbody.time << "\t " << nbody.arad << "\t " << nbody.steps << "\t " << 
	    total_energy << "\t " << -e_p << "\t " << e_k << "\t " << nbody.udda;
	  FileEnergy << "\t" << Nbody::omega(nbody.arad,nbody.omega_start,nbody.lambda_start);
	  FileEnergy << "\t" << Nbody::lambda(nbody.arad,nbody.omega_start,nbody.lambda_start);
	  FileEnergy << "\t" << Nbody::hubble(nbody.arad,nbody.omega_start,nbody.lambda_start) << endl;
	  FileMom << nbody.time << "\t " << nbody.arad << "\t " << nbody.steps << "\t " << 
	    sum_vx << "\t " << sum_vy << "\t " << sum_vz << "\t " <<
	    sum_px << "\t " << sum_py << "\t " << sum_pz << endl;
	}
    }
  else
    {
      double total_energy=0.5*(nbody.kinetic_energy_old+nbody.kinetic_energy)+
	nbody.potential_energy;
      if(nbody.steps >= 0)
      	{
	  FileEnergy << nbody.time <<  "\t " << nbody.steps << "\t " << 
	    total_energy << "\t " << -nbody.potential_energy << "\t " << 0.5*(nbody.kinetic_energy_old+nbody.kinetic_energy) << "\t " << endl;
	  FileMom << nbody.time << "\t " << nbody.arad << "\t " << nbody.steps << "\t " << 
	    sum_vx << "\t " << sum_vy << "\t " << sum_vz << "\t " <<
	    sum_px << "\t " << sum_py << "\t " << sum_pz << endl;
	}
    }
  nbody.kinetic_energy_old=nbody.kinetic_energy;
  nbody.potential_energy_old=nbody.potential_energy;
}
