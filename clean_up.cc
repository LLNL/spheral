#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void clean_up(Fractal_Memory& mem)
  {
    ofstream& FileFractal=mem.p_fractal->p_file->DUMPS;
    //--------------------------------------------------------------------------------------------------------------------------------
    // Delete all points, groups, chains and misc.
    //--------------------------------------------------------------------------------------------------------------------------------
    int counting=0;
    counting++;
    if(!mem.amnesia)
      return;
    //
    FileFractal << "number of everything entering clean_up "  << " " << Group::number_groups << " " << Point::number_points << "\n";
    int level_start=0;
    for(int level=level_start;level <= mem.level_max;level++)
      {
	for(auto p_group : mem.all_groups[level])
	  {
	    for(auto point_itr=p_group->list_new_points.begin();point_itr !=p_group->list_new_points.end();++point_itr)
	      {
		Point* p_whatever=*point_itr;
		delete [] p_whatever;
		p_whatever=0;
		counting++;
	      }
	    delete p_group;
	    counting++;
	    p_group=0;
	  }
      }
    Misc::DeleteMisc(); // Soliton
    mem.p_misc=0;
    counting++;
    mem.all_groups.clear();
    //
    FileFractal << "I have deleted this much " << counting << "\n";
    FileFractal << sizeof(int) << " int ";
    FileFractal << INT_MAX << " largest int ";
    FileFractal << sizeof(float) << " float ";
    FileFractal << FLT_MAX << " largest float " ;
    FileFractal << sizeof(double) << " double ";
    FileFractal << DBL_MAX << " largest double " << "\n";
    FileFractal << RAND_MAX << " rand_max " << "\n";
    FileFractal << sizeof(bool) << " bool ";
    FileFractal << sizeof(Group) << " group ";
    FileFractal << sizeof(Point) << " point ";
    FileFractal << sizeof(Fractal) << " fractal ";
    FileFractal << sizeof(Fractal_Memory) << " fractal memory ";
    FileFractal << sizeof(Misc) << " misc ";
    FileFractal << sizeof(Particle) << " Particle ";
    FileFractal << "\n";
    FileFractal << "number of everything exiting clean_up "  << " " << Group::number_groups << " " << Point::number_points << "\n";
  }
}
