#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void small_exceptions(Fractal_Memory& mem)
  {
    for(int level=1;level <= mem.global_level_max;level++)
      for(auto &pg : mem.all_groups[level])
	if(pg->list_points.size() == 27)
	  pg->set_buffer_group(false);
  }
}
