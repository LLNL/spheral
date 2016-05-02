#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void high_points(Group& group,Fractal& fractal,Misc& misc)
  {
    //    cout << "inside high_points " << " " << group.list_points.size() << endl;
    group.p_list_really_high.clear();
    int ni=0;
    for(vector <Point*>::const_iterator point_itr=group.list_points.begin();point_itr != group.list_points.end();++point_itr)
      {
	Point& point=**point_itr;
	if(check_high(point,fractal) && high_enough_level(point,group,fractal,misc))
	  {
	  point.set_it_is_high(true);
	  ni++;
	  }
      }
    //    cout << "ni= " << ni << endl;
    group.p_list_really_high.reserve(ni);
    for(vector <Point*>::const_iterator point_itr=group.list_points.begin();point_itr != group.list_points.end();++point_itr)
      {
	Point* p_point=*point_itr;
	if(p_point->get_it_is_high())
	  group.p_list_really_high.push_back(p_point);
      }
    group.set_number_high_points(ni);
  }
  bool high_enough_level(Point& point,Group& group,Fractal& fractal,Misc& misc)
  {
    vector <double> cen(3);
    vector <double> rad(3);
    bool square;
    if(fractal.get_number_masks() < 1) return true;
    const int lev=group.get_level();
    double move;
    if(fractal.get_padding() > 0 || misc.get_buffer_it())
      move=(1.0-pow(2.0,-lev))*misc.zoom;
    else
      move=0.0;
    if(move < 0.0) return 0;
    double scaling=(double)misc.grid_multiply;
    for(int m=0;m< fractal.get_number_masks();++m)
      {
	if(fractal.get_level_mask(m) > lev)
	  {
	    fractal.get_mask(m,cen,rad,square);
	    double d_x,d_y,d_z;
	    bool close=false;
	    double p_x=cen[0]*scaling;
	    double p_y=cen[1]*scaling;
	    double p_z=cen[2]*scaling;
	    double r_x=rad[0]*scaling;
	    double r_y=rad[1]*scaling;
	    double r_z=rad[2]*scaling;
	    if(min(min(r_x,r_y),r_z)-move < 0.0) return 0;
	    if(fractal.get_periodic())
	      {
		d_x=dist1((double)point.get_pos_point_x()-p_x,scaling);
		d_y=dist1((double)point.get_pos_point_y()-p_y,scaling);
		d_z=dist1((double)point.get_pos_point_z()-p_z,scaling);
	      }
	    else
	      {
		d_x=abs((double)point.get_pos_point_x()-p_x);
		d_y=abs((double)point.get_pos_point_y()-p_y);
		d_z=abs((double)point.get_pos_point_z()-p_z);
	      }
	    if(square)  
	      close= (d_x <= r_x-move) && (d_y <= r_y-move)  && (d_z <= r_z - move);
	    else
	      close=pow(d_x/(r_x-move),2)+pow(d_y/(r_y-move),2)+pow(d_z/(r_z-move),2) <= 1.0;
	    if(close) return true;
	  }
      }
    return false;
  }
  double dist1(const double& x,const double& y)
  {
    return min(min(abs(x),abs(x-y)),abs(x+y));
  }
}
