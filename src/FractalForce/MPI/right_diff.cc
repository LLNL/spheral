#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  bool right_diff(vector <int>& Va,vector <int>& Vb,vector <int>& VD)
  {
    bool ok=true;
    int ni=0;
    while(ok && ni < 3)
      {
	ok=Vb[ni]-Va[ni]==VD[ni];
	ni++;
      }
    return ok;
  }
}
