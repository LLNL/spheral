#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void update_rv(Fractal& fractal,const int& param,const double& const1,const double& const2)
  {
    ofstream& FileFractal=fractal.p_file->DUMPS;
    //    ofstream& FileFractal=fractal.p_file->FileFractal;
    vector <double> pos(3);
    vector <double> vel(3);
    vector <double> force(3);
    FileFractal << "update " << param << " " << const1 << " " << const2 << "\n";
    if(param == 0)
      {
	for(int n=0;n<fractal.get_number_particles();++n)
	  {
	    Particle* p=fractal.particle_list[n];
	    if(!p->get_real_particle())
	       continue;
	    p->get_pos(pos);
	    p->set_vel(pos);
	  }
      }
    else if(param == 1)
      {
	for(int n=0;n<fractal.get_number_particles();++n)
	  {
	    Particle* p=fractal.particle_list[n];
	    if(!p->get_real_particle())
	       continue;
	    p->get_pos(pos);
	    p->get_force(force);
	    pos[0]+=force[0]*const1;
	    pos[1]+=force[1]*const1;
	    pos[2]+=force[2]*const1;
	    p->set_pos(pos);
	  }
      }
    else if(param ==2)
      {
	for(int n=0;n<fractal.get_number_particles();++n)
	  {
	    Particle* p=fractal.particle_list[n];
	    if(!p->get_real_particle())
	       continue;
	    p->get_vel(vel);
	    p->get_force(force);
	    pos[0]=vel[0]+force[0]*const1;
	    pos[1]=vel[1]+force[1]*const1;
	    pos[2]=vel[2]+force[2]*const1;
	    vel[0]=force[0]*const2;
	    vel[1]=force[1]*const2;
	    vel[2]=force[2]*const2;
	    p->set_pos(pos);
	    p->set_vel(vel);
	  }
      }
    else if(param ==3)
      {
	double twopi=8.0*atan(1.0);
	double twopiacr=1.0/(twopi*const1);
	for(int n=0;n<fractal.get_number_particles();++n)
	  {
	    Particle* p=fractal.particle_list[n];
	    if(!p->get_real_particle())
	       continue;
	    p->get_vel(vel);
	    double dz=-twopiacr*sin(twopi*(vel[2]-const2));
	    pos[0]=vel[0];
	    pos[1]=vel[1];
	    pos[2]=vel[2]+dz;
	    vel[0]=0.0;
	    vel[1]=0.0;
	    vel[2]=dz;
	    p->set_pos(pos);
	    p->set_vel(vel);
	  }
      }
    else
      assert(0);
  }
}
