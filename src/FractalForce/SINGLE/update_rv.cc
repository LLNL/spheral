#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void update_rv(Fractal& fractal,const int& param,const double& const1,const double& const2)
  {
    vector <double> pos(3);
    vector <double> vel(3);
    vector <double> force(3);
    cout << "update " << param << " " << const1 << " " << const2 << endl;
    if(param == 0)
      {
	for(int n=0;n<fractal.get_number_particles();++n)
	  {
	    Particle* p=fractal.particle_list[n];
	    p->get_pos(pos);
	    p->set_vel(pos);
	  }
      }
    else if(param == 1)
      {
	for(int n=0;n<fractal.get_number_particles();++n)
	  {
	    Particle* p=fractal.particle_list[n];
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
    else
      assert(0);
  }
}
