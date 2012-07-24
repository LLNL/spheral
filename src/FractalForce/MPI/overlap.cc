#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  template <class T> bool overlap(vector <T>& xleft,vector <T>& xright,vector <T>& yleft,vector <T>& yright)
  {
    return (xleft[0] <= yright[0] && xright[0] >= yleft[0] && 
	    xleft[1] <= yright[1] && xright[1] >= yleft[1] &&
	    xleft[2] <= yright[2] && xright[2] >= yleft[2]);
  }
  template <class T> bool overlap(vector <T>& xleft,vector <T>& xright,vector <T>& yleftright)
  {
    return (xleft[0] <= yleftright[3] && xright[0] >= yleftright[0] && 
	    xleft[1] <= yleftright[4] && xright[1] >= yleftright[1] &&
	    xleft[2] <= yleftright[5] && xright[2] >= yleftright[2]);
  }
}
