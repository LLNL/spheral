#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void force_shear_at_point(Group& group, Fractal& fractal)
  {
    //    ofstream& FileFractal=fractal.p_file->DUMPS;
    //    ofstream& FileFractal=fractal.p_file->FileFractal;
    //    FileFractal << " enter force_shear_at_point " << "\n";
    const double conv=(double)(fractal.get_grid_length())*pow(2.0,group.get_level()-1);
    //
    for( vector <Point*>::const_iterator point_itr=group.list_points.begin();point_itr !=group.list_points.end();++point_itr)
      {
	Point& point=**point_itr;
// 	point.force_shear_point_make();
	const int rp=point.get_real_pointer();
	if(rp >= 0)
	  {
	    if(point.get_inside())
	      point.diff_force(conv);
	    else if(Point::corner[rp])
	      point.copy_force_shear_point_1();
	  }
      }
    for( vector <Point*>::const_iterator point_itr=group.list_points.begin();point_itr !=group.list_points.end();++point_itr)
      {
	Point& point=**point_itr;
	const int rp=point.get_real_pointer();
	if(rp >= 0)
	  {
	    if(!point.get_inside() && Point::edge[rp])
	      point.copy_force_shear_point_2(Point::cefc[rp]);
	  }
      }
    for( vector <Point*>::const_iterator point_itr=group.list_points.begin();point_itr !=group.list_points.end();++point_itr)
      {
	Point& point=**point_itr;
	const int rp=point.get_real_pointer();
	if(rp >=0)
	  {
	    if(!point.get_inside() && Point::face[rp])
	      point.copy_force_shear_point_4(Point::cefc[rp]);
	  }
      }
    //    FileFractal << " exit force shear at point " << "\n";
  }
}
