#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void potential_start(Group& group)
  {
    //---------------------------------------------------------------------------------------------------
    // We use potential from the mother group as initial conditions  and boundary conditions for the potential
    // We look at the 3x3x3 we made. First we do a straight copy at the 8 corner points. 
    // Then we interpolate at  the 12 edge points using the values at the 2 nearest corner.
    // Then we interpolate at the 6 face points using the 4 nearest edge points.
    // Finally we interpolate at the center using the 6 face points.
    //---------------------------------------------------------------------------------------------------
    for(auto &p_point : group.list_points)
      {
	if(Point::corner[p_point->get_real_pointer()])
	  p_point->copy_potential_point_1();
      }
    //
    for(auto &p_point : group.list_points)
      {
	int rp=p_point->get_real_pointer();
	if(Point::edge[rp])
	  p_point->copy_potential_point_2(Point::cefc[rp]);
      }
    //
    for(auto &p_point : group.list_points)
      {
	int rp=p_point->get_real_pointer();
	if(Point::face[rp])
	  p_point->copy_potential_point_4(Point::cefc[rp]);
      }
    //
    for(auto &p_point : group.list_points)
      {
	if(p_point->get_real_pointer() == 13)
	  p_point->copy_potential_point_6();
      }
  }
}
