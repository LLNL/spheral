#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void potential_start(Group& group)
  {
    //    cerr << "enter pot start " << &group << "\n";
    //---------------------------------------------------------------------------------------------------
    // We use potential from the mother group as initial conditions  and boundary conditions for the potential
    // We look at the 3x3x3 we made. First we do a straight copy at the 8 corner points. 
    // Then we interpolate at  the 12 edge points using the values at the 2 nearest corner.
    // Then we interpolate at the 6 face points using the 4 nearest edge points.
    // Finally we interpolate at the center using the 6 face points.
    //---------------------------------------------------------------------------------------------------
    for(vector <Point*>::const_iterator point_itr=group.list_points.begin();point_itr != group.list_points.end();++point_itr)
      {
	Point* p_point=*point_itr;
	if(Point::corner[p_point->get_real_pointer()])
	  p_point->copy_potential_point_1();
      }
    //
    for(vector <Point*>::const_iterator point_itr=group.list_points.begin();point_itr != group.list_points.end();++point_itr)
      {
	Point* p_point=*point_itr;
	int rp=p_point->get_real_pointer();
	if(Point::edge[rp])
	  p_point->copy_potential_point_2(Point::cefc[rp]);
      }
    //
    for(vector <Point*>::const_iterator point_itr=group.list_points.begin();point_itr != group.list_points.end();++point_itr)
      {
	Point* p_point=*point_itr;
	int rp=p_point->get_real_pointer();
	if(Point::face[rp])
	  p_point->copy_potential_point_4(Point::cefc[rp]);
      }
    //
    for(vector <Point*>::const_iterator point_itr=group.list_points.begin();point_itr != group.list_points.end();++point_itr)
      {
	Point* p_point=*point_itr;
	if(p_point->get_real_pointer() == 13)
	  p_point->copy_potential_point_6();
      }
    //
    //    cerr << "exit pot start" << "\n";
  }
}
