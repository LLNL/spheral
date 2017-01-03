#include "libs.hh"
#include "classes.hh"
#include "nbody.hh"
#include "headers.hh"
void take_a_step(Nbody& nbody, Fractal& fractal)
{
  if(nbody.periodic)
    {
      double parad=pow(nbody.arad,nbody.pexp);
      double dpda=nbody.pexp*parad/nbody.arad;
      double parad_half=parad+0.5*nbody.step_length;
      double arad_half=pow(parad_half,1.0/nbody.pexp);
      double dadt=arad_half*Nbody::hubble(arad_half,nbody.omega_start,nbody.lambda_start);
      double v_const=1.0-2.0*nbody.step_length/nbody.arad/dpda;
      double f_const=nbody.step_length/pow(nbody.arad,4)/Nbody::hubble(arad_half,nbody.omega_start,nbody.lambda_start)/dpda;
      double dt=nbody.step_length/dadt/dpda;
      for(int n=0;n < nbody.number_particles;++n)
	{
	  nbody.vel_x[n]=nbody.vel_x[n]*v_const+fractal.get_force_x(n)*f_const;
	  nbody.vel_y[n]=nbody.vel_y[n]*v_const+fractal.get_force_y(n)*f_const;
	  nbody.vel_z[n]=nbody.vel_z[n]*v_const+fractal.get_force_z(n)*f_const;
	  nbody.pos_x[n]=fractal.get_pos_x(n)+nbody.vel_x[n]*dt;
	  nbody.pos_y[n]=fractal.get_pos_y(n)+nbody.vel_y[n]*dt;
	  nbody.pos_z[n]=fractal.get_pos_z(n)+nbody.vel_z[n]*dt;
	  //	 nbody.highest_level[n]=fractal.get_highest_level(n);
	}
      nbody.time+=nbody.step_length/dadt/dpda;
      if(!fractal.get_parameters_bool(513))      // not startup
	parad+=nbody.step_length;
      nbody.arad=pow(parad,1.0/nbody.pexp);
    }
  else
    {
      double dt=nbody.step_length;
      for(int n=0;n < nbody.number_particles;++n)
	{
	  nbody.vel_x[n]+=fractal.get_force_x(n)*dt;
	  nbody.vel_y[n]+=fractal.get_force_y(n)*dt;
	  nbody.vel_z[n]+=fractal.get_force_z(n)*dt;
	  nbody.pos_x[n]+=nbody.vel_x[n]*dt;
	  nbody.pos_y[n]+=nbody.vel_y[n]*dt;
	  nbody.pos_z[n]+=nbody.vel_z[n]*dt;
	  //	 nbody.highest_level[n]=fractal.get_highest_level(n);
	}
    }
}
