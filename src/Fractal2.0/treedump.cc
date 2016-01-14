#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void tree_dump(Fractal_Memory& FM)
  {
    vector <int>pos(3);
    for(int level=0;level <= FM.level_max;level++)
      {
	int gnumber=0;
	for(vector <Group*>::const_iterator group_itr=fractal_memory.all_groups[level].begin();
	    group_itr!=fractal_memory.all_groups[level].end();group_itr++)
	  {
	    int pnumber=0;
	    Group& group=**group_itr;
	    for( vector <Point*>::const_iterator point_itr=group.list_points.begin();point_itr !=group.list_points.end();++point_itr)
	      {
		Point& point=**point_itr;
		int rp=point.get_real_pointer();
		point.get_pos_point(pos);
		int npart=point.list_particles.size();
		FM.p_file->FileDump << " S" << FM.steps << " G" << gnumber << " P" << rp << " " << npart;
		FM.p_file->FileDump << " " << pos[0]  << " " << pos[1]  << " " << pos[2] << "\n";
		pnumber++;
	      }
	    gnumber++;
	  }
      }
  }
}
