#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void force_at_point(Group& group, Fractal& fractal)
  {
    // worry about group 1 at the edge for isolated BC.
    //    ofstream& FileFractal=fractal.p_file->DUMPS;
    //    ofstream& FileFractal=fractal.p_file->FileFractal;
    //    if(fractal.get_debug()) FileFractal << " enter force at point " << &group << " " << group.get_level() << "\n";
    const double conv=(double)(fractal.get_grid_length())*pow(2.0,group.get_level()-1);
    //
    if(group.get_level() == 0)
      {
	for( vector <Point*>::const_iterator point_itr=group.list_points.begin();point_itr !=group.list_points.end();++point_itr)
	  {
	    Point& point=**point_itr;
	    const int rp=point.get_real_pointer();
	    if(rp < 0)
	      continue;
	    if(point.get_inside())
	      point.diff_pot(conv);
	    else
	      point.diff_pot_careful(conv);
	    //	    point.dumppf();
	  }
	return;
      }
    for( vector <Point*>::const_iterator point_itr=group.list_points.begin();point_itr !=group.list_points.end();++point_itr)
      {
	Point& point=**point_itr;
	const int rp=point.get_real_pointer();
	if(rp < 0)
	  continue;
	if(point.get_inside())
	  point.diff_pot(conv);
	else if(Point::corner[rp])
	  point.copy_force_point_1();
      }
    //    assert(0);
    for( vector <Point*>::const_iterator point_itr=group.list_points.begin();point_itr !=group.list_points.end();++point_itr)
      {
	Point& point=**point_itr;
	const int rp=point.get_real_pointer();
	if(rp < 0)
	  continue;
	if(!point.get_inside() && Point::edge[rp])
	  point.copy_force_point_2(Point::cefc[rp]);
      }
    for( vector <Point*>::const_iterator point_itr=group.list_points.begin();point_itr !=group.list_points.end();++point_itr)
      {
	Point& point=**point_itr;
	const int rp=point.get_real_pointer();
	if(rp < 0)
	  continue;
	if(!point.get_inside() && Point::face[rp])
	  point.copy_force_point_4(Point::cefc[rp]);
      }
    /*
    for( vector <Point*>::const_iterator point_itr=group.list_points.begin();point_itr !=group.list_points.end();++point_itr)
      {
	Point& point=**point_itr;
	point.dumppf();
      }
    */
    //    if(fractal.get_debug()) FileFractal << " exit force at point " << "\n";
  }
}
