#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void hypre_struct_compact_boxes(vector <Point*>& p_points_left,vector <Point*>& p_points_right,Fractal& fractal,Fractal_Memory& mem,const int& level)
  {
    const unsigned int zoom=Misc::pow(2,fractal.get_level_max()-level);
    const int boxes=p_points_left.size();
    int nj=0;
    int nk=1;
    int nb=0;
    vector <int>posjl(3);
    vector <int>posjr(3);
    vector <int>poskl(3);
    vector <int>poskr(3);
    for(int ni=0;ni<boxes-1;ni++)
      {
	p_points_left[nj]->get_pos_point(posjl);
	p_points_right[nj]->get_pos_point(posjr);
	p_points_left[nk]->get_pos_point(poskl);
	p_points_right[nk]->get_pos_point(poskr);
	bool isit = is_this_a_box(posjl,posjr,poskl,poskr,zoom);
	if(isit)
	  {

	  }
	else
	  {
	    nj++;
	  }
      }
  }
  bool is_this_a_box(vector <int>& posjl,vector <int>& posjr,vector <int>& poskl,vector <int>& poskr,const int& zoom)
  {
    bool equal_x=poskl[0]==posjl[0];
    bool equal_y=poskl[1]==posjl[1];
    bool equal_z=poskl[2]==posjl[2];
    if(equal_x && equal_y)
      {
	return poskr[0]==posjr[0] && poskr[1]==posjr[1] && 
	  (poskl[2]-posjr[2]==zoom || posjl[2]-poskr[2]==zoom);
      }
    else if(equal_x && equal_z)
      {
	return poskr[0]==posjr[0] && poskr[2]==posjr[2] &&
	  (poskl[1]-posjr[1]==zoom || posjl[1]-poskr[1]==zoom);
      }
    else if(equal_y && equal_z)
      {
	return poskr[1]==posjr[1] && poskr[2]==posjr[2] &&
	  (poskl[0]-posjr[0]==zoom || posjl[0]-poskr[0]==zoom);
      }
    return false;
  }
}
