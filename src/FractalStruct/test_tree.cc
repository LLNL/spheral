#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  bool test_tree(Fractal_Memory& fractal_memory,Fractal& fractal)
  {
    fractal.timing(-1,16);
    ofstream& FileFractal=fractal.p_file->DUMPS;
    bool badd=false;
    for(int level=1;level <= fractal.get_level_max();level++)
      {
	for(auto &pgroup : fractal_memory.all_groups[level])
	  {
	    Group& group=*pgroup;
	    sort3_list(*pgroup,0); // sort the points in x,y,z and test for dupes
	    if(test_group(*pgroup)) FileFractal << "badd group " << &group << group.get_level() << "\n";
	  }
      }
    fractal.timing(1,16);
    return badd;
  }
}
