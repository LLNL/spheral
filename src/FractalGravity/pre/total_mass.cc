#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  double total_mass(Fractal& frac)
  {
    double sum=0.0;
    for(int ni=0;ni<frac.get_number_particles();++ni)
      sum+=frac.particle_list[ni]->get_mass();
    return sum;
  }
}
