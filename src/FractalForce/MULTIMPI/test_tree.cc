#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  bool test_tree(Fractal_Memory& fractal_memory,Fractal& fractal)
  {
    ofstream& FileFractal=fractal.p_file->FileFractal;
    bool badd=false;
    for(int level=1;level <= fractal.get_level_max();level++)
      {
	for(vector <Group*>::const_iterator group_itr=fractal_memory.all_groups[level].begin();
	    group_itr!=fractal_memory.all_groups[level].end();group_itr++)
	  {
	    Group& group=**group_itr;
	    sort3_list(group,0); // sort the points in x,y,z and test for dupes
	    if(test_group(group)) FileFractal << "badd group " << &group << group.get_level() << endl;
	  }
      }
    return badd;
  }
}
