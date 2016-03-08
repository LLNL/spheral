#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void density_edge(Group& group, Misc& misc)
  {
    for(list <Point*>::const_iterator point_itr=group.list_points.begin();point_itr != group.list_points.end();++point_itr)
      {
	Point& point=**point_itr;
	int rp=point.get_real_pointer();
	if(Misc::is_it_corner(rp))
	  {
	    Point* p_p=point.get_point_pointer();
	    point.set_density_point(p_p->get_density_point());
	  }
      }
    //
    for(list<Point*>::const_iterator point_itr=group.list_points.begin();point_itr !=group.list_points.end();++point_itr)
      {
	Point& point=**point_itr;
	int rp=point.get_real_pointer();
	if(Misc::is_it_edge(rp))
	  {
	    int rp_x=rp % 3;
	    int rp_y=rp/3 % 3;
	    Point* p_up=0;
	    Point* p_down=0;
	    //
	    if(rp_x % 2 == 1)
	      {
		p_up=point.get_point_up_x_0();
		p_down=point.get_point_down_x_0();
	      }
	    else if(rp_y % 2 == 1)
	      {
		p_up=point.get_point_up_y_0();
		p_down=point.get_point_down_y_0();
	      }
	    else
	      {
		p_up=point.get_point_up_z_0();
		p_down=point.get_point_down_z_0();
	      }
	    point.copy_density_point(*p_up,*p_down);
	  }
      }
    for(list<Point*>::const_iterator point_itr=group.list_points.begin();point_itr !=group.list_points.end();++point_itr)
      {
	Point& point=**point_itr;
	int rp=point.get_real_pointer();
	if(Misc::is_it_face(rp))
	  {
	    int rp_y=rp/3 % 3;
	    int rp_z=rp/9;
	    Point* p_up_0=0;
	    Point* p_down_0=0;
	    Point* p_up_1=0;
	    Point* p_down_1=0;
	    if(rp_z % 2 == 0)
	      {
		p_up_0=point.get_point_up_x_0();
		p_up_1=point.get_point_up_y_0();
		p_down_0=point.get_point_down_x_0();
		p_down_1=point.get_point_down_y_0();
	      }
	    else if(rp_y % 2 == 0)
	      {
		p_up_0=point.get_point_up_x_0();
		p_up_1=point.get_point_up_z_0();
		p_down_0=point.get_point_down_x_0();
		p_down_1=point.get_point_down_z_0();
	      }
	    else
	      {
		p_up_0=point.get_point_up_y_0();
		p_up_1=point.get_point_up_z_0();
		p_down_0=point.get_point_down_y_0();
		p_down_1=point.get_point_down_z_0();
	      }
	    point.copy_density_point(*p_up_0,*p_down_0,*p_up_1,*p_down_1);
	  }
      }
    // no need to do rp=13 points, they are always inside
    for(list <Point*>::const_iterator point_itr=group.list_points.begin();point_itr != group.list_points.end();++point_itr)
      {
	if((*point_itr)->get_inside()) (*point_itr)->set_density_point(0.0);
      }
  }
}
