#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  int findPointByPos(vector <Point*>& plist,Point* psend,FILE* FF)
  {
    vector <Point*>::iterator found=std::lower_bound(plist.begin(),plist.end(),psend,LesserPoint);
    if(found == plist.end())
      {
	fprintf(FF," not found \n");
	return -1;
      }
    int label= (found-plist.begin());
    return label;
  }
}
