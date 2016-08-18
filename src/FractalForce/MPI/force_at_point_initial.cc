#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void force_at_point_initial(Group& group, Fractal& fractal)
  {
    using namespace std;
    ofstream& FileFractal=fractal.p_file->FileFractal;
    // worry about group 1 at the edge for isolated BC.
    FileFractal << " enter force at point initial" << endl;
    //
    for( vector <Point*>::const_iterator point_itr=group.list_points.begin();point_itr !=group.list_points.end();++point_itr)
      {
	Point& point=**point_itr;
	const int rp=point.get_real_pointer();
	if(rp >= 0)
	  {
	    if(Point::corner[rp])
	      point.copy_force_point_1();
	  }
      }
    for( vector <Point*>::const_iterator point_itr=group.list_points.begin();point_itr !=group.list_points.end();++point_itr)
      {
	Point& point=**point_itr;
	const int rp=point.get_real_pointer();
	if(rp >= 0)
	  {
	    if(Point::edge[rp])
	      point.copy_force_point_2(Point::cefc[rp]);
	  }
      }
    for( vector <Point*>::const_iterator point_itr=group.list_points.begin();point_itr !=group.list_points.end();++point_itr)
      {
	Point& point=**point_itr;
	const int rp=point.get_real_pointer();
	if(rp >=0)
	  {
	    if(Point::face[rp])
	      point.copy_force_point_4(Point::cefc[rp]);
	  }
      }
    for( vector <Point*>::const_iterator point_itr=group.list_points.begin();point_itr !=group.list_points.end();++point_itr)
      {
	Point& point=**point_itr;
	const int rp=point.get_real_pointer();
	if(rp == 13)
	  point.copy_force_point_6();
      }
    FileFractal << " exit force at point " << endl;
  }
}
