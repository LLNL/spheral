#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void boxes(Fractal& fractal,Group& group,vector < vector <int> >& boxab)
  {
    boxab.clear();
    int level=group.get_level();
    int delta=Misc::pow(2,fractal.get_level_max()-level);
    vector <int> posab(6);
    vector <int> posb(3);
    for(vector<Point*>::const_iterator point_itr=group.list_points.begin();point_itr !=group.list_points.end();++point_itr)
      {
	Point& point=**point_itr;
	if(point.get_real_pointer() != 13)
	  continue;
	point.get_pos_point(posb);
	posab[0]=posb[0]/delta-1;
	posab[1]=posb[1]/delta-1;
	posab[2]=posb[2]/delta-1;
	posab[3]=posab[0]+2;
	posab[4]=posab[1]+2;
	posab[5]=posab[2]+2;
	boxab.push_back(posab);
	/*
	Point* up_x=point.get_point_up_x();
	Point* up_y=point.get_point_up_y();
	Point* up_z=point.get_point_up_z();
	*/
      }
  }
}
