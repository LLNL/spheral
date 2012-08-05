#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void force_at_point(Group& group, Fractal& fractal)
  {
    //  using namespace std;
    // worry about group 1 at the edge for isolated BC.
    //    if(fractal.get_debug()) cout << " enter force at point " << &group << " " << group.get_level() << endl;
    const double conv=(double)(fractal.get_grid_length())*pow(2.0,group.get_level()-1);
    //
    for( vector <Point*>::const_iterator point_itr=group.list_points.begin();point_itr !=group.list_points.end();++point_itr)
      {
	Point& point=**point_itr;
	const int rp=point.get_real_pointer();
	if(rp >= 0)
	  {
	    if(point.get_inside())
	      point.diff_pot(conv);
	    else if(Point::corner[rp])
	      point.copy_force_point_1();
	  }
      }
    for( vector <Point*>::const_iterator point_itr=group.list_points.begin();point_itr !=group.list_points.end();++point_itr)
      {
	Point& point=**point_itr;
	const int rp=point.get_real_pointer();
	if(rp >= 0)
	  {
	    if(!point.get_inside() && Point::edge[rp])
	      {
		point.copy_force_point_2(Point::cefc[rp]);
	      }
	  }
      }
    for( vector <Point*>::const_iterator point_itr=group.list_points.begin();point_itr !=group.list_points.end();++point_itr)
      {
	Point& point=**point_itr;
	const int rp=point.get_real_pointer();
	if(rp >=0)
	  {
	    if(!point.get_inside() && Point::face[rp])
	      {
		point.copy_force_point_4(Point::cefc[rp]);
	      }
	  }
      }
    //    if(fractal.get_debug()) cout << " exit force at point " << endl;
  }
}
