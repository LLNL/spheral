#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  int which_element(vector <Point*>& vec,const int& x,const int& y,const int& z,ofstream& FF)
  {
    vector <int>pos(3);
    int length=vec.size();
    for(int ni=0;ni<length;ni++)
      {
	vec[ni]->get_pos_point(pos);
	if(pos[0] == x && pos[1] == y && pos[2] == z)
	  return ni;
      }
    FF << " did not find it " << x << " " << y << " " << z << endl;
    for(int ni=0;ni<length;ni++)
      {
	vec[ni]->get_pos_point(pos);
	FF << ni << " " << pos[0] << " " << pos[1] << " " << pos[2] << endl;
      }
    //    assert(0);
    return -1;
  }
}
