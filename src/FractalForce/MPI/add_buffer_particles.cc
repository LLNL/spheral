#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void add_buffer_particles(Fractal_Memory& mem,Fractal& frac)
  {
    if(!mem.periodic)
      return;
    int length=mem.grid_length;
    double Rlow=-2.0/static_cast<double>(length);
    double Rhigh=1.0-Rlow;
    vector <double>pos(3);
    for(int particle=0; particle < frac.get_number_particles(); ++particle)
      {
	Particle* P=frac.particle_list[particle];
	P->get_pos(pos);
	for(int nz=-1;nz<=1;nz++)
	  {
	    double posz=pos[2]+nz;
	    int okz=posz > Rlow && posz < Rhigh;
	    if(!okz)
	      continue;
	    for(int ny=-1;ny<=1;ny++)
	      {
		double posy=pos[1]+ny;
		int oky=posy > Rlow && posy < Rhigh;
		if(!oky)
		  continue;
		for(int nx=-1;nx<=1;nx++)
		  {
		    double posx=pos[0]+nx;
		    int okx=posx > Rlow && posx < Rhigh;
		    if(!okx || (nz==0 && ny==0 && nx==0 ))
		      continue;
		    double m=P->get_mass();
		    Particle* Pb=new Particle;
		    Pb->set_pos(posx,posy,posz);
		    Pb->set_mass(m);
		    frac.particle_list.push_back(Pb);
		  }
	      }
	  }
      }
  }
}
