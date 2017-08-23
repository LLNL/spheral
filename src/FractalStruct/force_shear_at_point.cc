#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void force_shear_at_point(Group& group, Fractal& fractal)
  {
    const double conv=(double)(fractal.get_grid_length())*pow(2.0,group.get_level()-1);
    for(auto &p : group.list_points)
      {
	const int rp=p->get_real_pointer();
	if(rp < 0)
	  continue;
	if(p->get_inside())
	  p->diff_force(conv);
	else if(Point::corner[rp])
	  p->copy_force_shear_point_1();
      }
    for(auto &p : group.list_points)
      {
	const int rp=p->get_real_pointer();
	if(rp >= 0)
	  if(!p->get_inside() && Point::edge[rp])
	    p->copy_force_shear_point_2(Point::cefc[rp]);
      }
    for(auto &p : group.list_points)
      {
	const int rp=p->get_real_pointer();
	if(rp >=0)
	  if(!p->get_inside() && Point::face[rp])
	    p->copy_force_shear_point_4(Point::cefc[rp]);
      }
  }
}
