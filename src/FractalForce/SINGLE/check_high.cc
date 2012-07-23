#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  bool check_high(Point& point,Fractal& fractal)
  {
    //--------------------------------------------------------------------------------------------------------------------------------
    // Guess what, as Nina says.
    //--------------------------------------------------------------------------------------------------------------------------------
    if(point.list_particles.size() < fractal.get_minimum_number()) return false;
    unsigned int np=0;
    for(vector<Particle*>::const_iterator particle_itr=point.list_particles.begin();particle_itr !=point.list_particles.end();++particle_itr)
      {
	if((*particle_itr)->get_p_highest_level_group() != 0) np++;
	if(np  >= fractal.get_minimum_number()) return true;
      }
    return false;
  }
}
