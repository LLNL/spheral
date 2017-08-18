#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void force_at_point(Group& group, Fractal& fractal)
  {
    fractal.timing(-1,7);
    const double conv=(double)(fractal.get_grid_length())*pow(2.0,group.get_level()-1);
    if(group.get_level() == 0)
      {
	for(Point* &p_point : group.list_points)
	  {
	    const int rp=p_point->get_real_pointer();
	    if(rp < 0)
	      continue;
	    if(p_point->get_inside())
	      p_point->diff_pot(conv);
	    else
	      p_point->diff_pot_careful(conv);
	    //	    p_point->dumppf();
	  }
	return;
      }
    for(Point* &p_point : group.list_points)
      {
	const int rp=p_point->get_real_pointer();
	if(rp < 0)
	  continue;
	if(p_point->get_inside() && !p_point->get_trouble())
	  p_point->diff_pot(conv);
	else if(Point::corner[rp])
	  p_point->copy_force_point_1();
      }
    for(Point* &p_point : group.list_points)
      {
	const int rp=p_point->get_real_pointer();
	if(rp < 0)
	  continue;
	if(!p_point->get_inside() && Point::edge[rp])
	  p_point->copy_force_point_2(Point::cefc[rp]);
      }
    for(Point* &p_point : group.list_points)
      {
	const int rp=p_point->get_real_pointer();
	if(rp < 0)
	  continue;
	if(!p_point->get_inside() && Point::face[rp])
	  p_point->copy_force_point_4(Point::cefc[rp]);
      }
    fractal.timing(1,7);
  }
}
