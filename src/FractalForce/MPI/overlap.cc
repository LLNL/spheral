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
}
namespace FractalSpace
{
  template bool overlap(vector <int>& xleft,vector <int>& xright,vector <int>& yleft,vector <int>& yright);
  template bool overlap(vector <double>& xleft,vector <double>& xright,vector <double>& yleft,vector <double>& yright);
}
namespace FractalSpace
{
  template <class T> bool overlap(vector <T>& xleft,vector <T>& xright,vector <T>& box)
  {
    return (xleft[0] <= box[1] && xright[0] >= box[0] && 
	    xleft[1] <= box[3] && xright[1] >= box[2] &&
	    xleft[2] <= box[5] && xright[2] >= box[4]);
  }
}
namespace FractalSpace
{
  template bool overlap(vector <int>& xleft,vector <int>& xright,vector <int>& box);
  template bool overlap(vector <double>& xleft,vector <double>& xright,vector <double>& box);
}
namespace FractalSpace
{
  template <class T> bool overlap(vector <T>& xvec,vector <T>& box)
  {
    return (xvec[0] <= box[1] && xvec[0] >= box[0] && 
	    xvec[1] <= box[3] && xvec[1] >= box[2] &&
	    xvec[2] <= box[5] && xvec[2] >= box[4]);
  }
}
namespace FractalSpace
{
  template bool overlap(vector <int>& xvec,vector <int>& box);
  template bool overlap(vector <double>& xvec,vector <double>& box);
}
