#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void check_for_edge_trouble(Fractal& fractal)
  { 
    ofstream& FileFractal=fractal.p_file->DUMPS;
    //    ofstream& FileFractal=fractal.p_file->FileFractal;
    fractal.timing(-1,2);
    //--------------------------------------------------------------------------------------------------------------------------------
    // Round off errors can cause trouble at the edge, move points a little
    //--------------------------------------------------------------------------------------------------------------------------------
    FileFractal << "edge trouble " << "\n";
    double eps=DBL_EPSILON;
    vector <double>pos(3);
    int outsiders=0;
    for(int part=0; part < fractal.get_number_particles();part++)
      {
	Particle* p=fractal.particle_list[part];
	if(p->get_p_highest_level_group() == 0)
	  continue;
	p->get_pos(pos);
	bool outside=pos[0] >= 1.0 || pos[0] <=0.0 ||
	  pos[1] >= 1.0 || pos[1] <=0.0 ||
	  pos[2] >= 1.0 || pos[2] <=0.0;
	if(!outside) continue;
	outsiders++;
	if(pos[0] >= 1.0)
	  pos[0]-=eps;
	else if(pos[0] <= 0.0) 
	  pos[0]+=eps;
	if(pos[1] >= 1.0) 
	  pos[1]-=eps;
	else if(pos[1] <= 0.0) 
	  pos[1]+=eps;
	if(pos[2] >= 1.0) 
	  pos[2]-=eps;
	else if(pos[2] <= 0.0) 
	  pos[2]+=eps;
	p->set_pos(pos);
      }
    FileFractal << " Total Outsiders " << outsiders << "\n";
    fractal.timing(1,2);
  }
}
