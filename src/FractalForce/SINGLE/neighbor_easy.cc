#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void neighbor_easy(vector <Point*>& p)
  {
    for( int rc=0;rc<3;rc++)
      {
	for(int rb=0;rb<3;rb++)
	  {
	    for(int ra=0;ra<3;ra++)
	      {
		int r1=ra+3*rb+9*rc;
		if(p[r1] != 0)
		  {
		    if(ra <2)
		      p[r1]->set_point_up_x(p[r1+1]);
		    if(rb <2)
		      p[r1]->set_point_up_y(p[r1+3]);
		    if(rc <2)
		      p[r1]->set_point_up_z(p[r1+9]);
		  }
	      }
	  }
      }
    for( int rc=0;rc<3;rc++)
      {
	for(int rb=0;rb<3;rb++)
	  {
	    for(int ra=0;ra<3;ra++)
	      {
		int r1=ra+3*rb+9*rc;
		if(p[r1] != 0)
		  p[r1]->down_from_up();
	      }
	  }
      }
  }
}
