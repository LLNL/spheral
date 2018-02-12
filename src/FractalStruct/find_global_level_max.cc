#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void find_global_level_max(Fractal_Memory& mem)
  {
    int highest_level_global=0;
    for(int level=0;level <= mem.level_max;level++)
      if(mem.all_groups[level].size() > 0)
	highest_level_global=level;
    vector <int>highest(1,highest_level_global);
    mem.p_mess->Find_Max_INT(highest,1);
    mem.global_level_max=highest[0];
  }
}
