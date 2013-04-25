#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  int find_global_level_max(Fractal_Memory& mem,Fractal& frac)
  {
    int highest_level_global=0;
    for(int level=0;level <= frac.get_level_max();level++)
      {
	if(mem.all_groups[level].size() > 0)
	  highest_level_global=level;
      }
    if(!mem.MPIrun)
      return highest_level_global;
    frac.timing(-1,41);
    mem.p_mess->Full_Stop();
    frac.timing(1,41);
    int* highest= new int[1];
    highest[0]=highest_level_global;
    mem.p_mess->Find_Max_INT(highest,1);
    highest_level_global=highest[0];
    delete [] highest;
    return highest_level_global;
  }
}
