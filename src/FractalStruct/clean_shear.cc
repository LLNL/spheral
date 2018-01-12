#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
//
namespace FractalSpace
{
  void clean_shear(Fractal_Memory& mem)
  {
    for(int level=0;level<=mem.p_fractal->get_level_max();level++)
      for(auto &pg : mem.all_groups[level])
	for(auto &p : pg->list_points)
	  p->clean_shear();
  }
}
