#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void add_pseudo_particles(Fractal_Memory& mem,Fractal& frac)
  {
    frac.p_file->FileFractal << " entered into add pseudo particles " << endl;
    if(!mem.periodic)
      return;
    ofstream& FF=mem.p_file->FileFractal;
    int length=mem.grid_length;
    double Rdelta=1.0/static_cast<double>(length);
    double Rlow=-2.0*Rdelta;
    double Rhigh=1.0+Rdelta;
    vector <double>pos(3);
    for(int particle=0; particle < frac.get_number_particles(); ++particle)
      {
	Particle* P=frac.particle_list[particle];
	//	FF << " add " << particle << " " ;
	//	P->dump(mem.p_file->FileFractal);
	P->get_pos(pos);
	for(int nz=-1;nz<=1;nz++)
	  {
	    double posz=pos[2]+nz;
	    bool okz=posz > Rlow && posz < Rhigh;
	    if(!okz)
	      continue;
	    for(int ny=-1;ny<=1;ny++)
	      {
		double posy=pos[1]+ny;
		bool oky=posy > Rlow && posy < Rhigh;
		if(!oky)
		  continue;
		for(int nx=-1;nx<=1;nx++)
		  {
		    double posx=pos[0]+nx;
		    bool okx=posx > Rlow && posx < Rhigh;
		    if(!okx || (nz==0 && ny==0 && nx==0 ))
		      continue;
		    double m=P->get_mass();
		    Particle* Pb=new Particle;
		    Pb->set_pos(posx,posy,posz);
		    Pb->set_mass(m);
		    Pb->set_real_particle(false);
		    frac.particle_list.push_back(Pb);
		    //		    Pb->dump(mem.p_file->FileFractal);
		  }
	      }
	  }
      }
    int parts=frac.particle_list.size();
    frac.set_number_particles(parts);
  }
}
