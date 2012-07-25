#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void left_right(Fractal& frac,vector <double>& pos_left,vector <double>& pos_right)
  {
    pos_left.assign(3,1.0e30);
    pos_right.assign(3,-1.0e30);
    vector <double> pos(3);
    for(int particle=0; particle < frac.get_number_particles(); ++particle)
      {
	frac.particle_list[particle]->get_pos(pos);
	for(int ni=0;ni<3;ni++)
	  {
	    pos_left[ni]=min(pos_left[ni],pos[ni]);
	    pos_right[ni]=max(pos_right[ni],pos[ni]);
	  }
      }
  }	
}
