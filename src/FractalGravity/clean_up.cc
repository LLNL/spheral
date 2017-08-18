#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void clean_up(Fractal_Memory& mem,Misc& misc,Fractal& fractal_ghost)
  {
    //--------------------------------------------------------------------------------------------------------------------------------
    // Delete all points, groups, chains and misc.
    //--------------------------------------------------------------------------------------------------------------------------------
    int counting=0;
    for(int particle=0;particle < fractal_ghost.get_number_particles();++particle)
      {
      delete fractal_ghost.particle_list[particle];
      counting++;
      }
    delete &fractal_ghost;
    counting++;
    if(!mem.amnesia)
      return;
    //
    cout << "number of everything entering clean_up "  << " " << Group::number_groups << " " << Point::number_points << endl;
    int level_start=0;
    if(mem.remember_points)
      level_start=1;
    for(int level=level_start;level <= mem.level_max;level++)
      {
	for(vector <Group*>::const_iterator group_itr=mem.all_groups[level].begin();
	    group_itr!=mem.all_groups[level].end();group_itr++)
	  {
	    Group* p_group=*group_itr;
	    for(vector <Point*>::const_iterator point_itr=p_group->list_new_points.begin();point_itr !=p_group->list_new_points.end();++point_itr)
	      {
		Point* p_whatever=*point_itr;
		delete [] p_whatever;
		p_whatever=0;
		counting++;
	      }
	    delete  p_group;
	    counting++;
	    p_group=0;
	  }
      }
    if(mem.remember_points)
      {
	Group* p_group_0=misc.p_group_0;
	for(unsigned int ni=0;ni<(p_group_0->list_points).size();ni++)
	  p_group_0->list_points[ni]->clean_point();
      }
    delete  &misc;
    mem.p_misc=0;
    counting++;
    mem.all_groups.clear();
    //
    cout << "I have deleted this much " << counting << endl;
    cout << sizeof(int) << " int ";
    cout << INT_MAX << " largest int ";
    cout << sizeof(float) << " float ";
    cout << FLT_MAX << " largest float " ;
    cout << sizeof(double) << " double ";
    cout << DBL_MAX << " largest double " << endl;
    cout << RAND_MAX << " rand_max " << endl;
    cout << sizeof(bool) << " bool ";
    cout << sizeof(Group) << " group ";
    cout << sizeof(Point) << " point ";
    cout << sizeof(Fractal) << " fractal ";
    cout << sizeof(Fractal_Memory) << " fractal memory ";
    cout << sizeof(Misc) << " misc ";
    cout << sizeof(Particle) << " Particle ";
    cout << endl;
    cout << "number of everything exiting clean_up "  << " " << Group::number_groups << " " << Point::number_points << endl;
  }
}
