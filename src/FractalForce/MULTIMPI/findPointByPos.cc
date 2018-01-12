#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  int findPointByPos(vector <Point*>& plist,Point* psend,FILE* FF)
  {
    //    FF << "find this " << psend->get_pos_point(0) << " " << psend->get_pos_point(1) << " " << psend->get_pos_point(2) << "\n";
    vector <Point*>::iterator found=std::lower_bound(plist.begin(),plist.end(),psend,LesserPoint);
    if(found == plist.end())
      {
	fprintf(FF," not found \n");
	return -1;
      }
    int label= (found-plist.begin());
    //    FF << "found this " << label << " " << plist[label] << " " << plist[label]->get_pos_point(0) << " " << plist[label]->get_pos_point(1) << " " << plist[label]->get_pos_point(2) << "\n";
    return label;
  }
}
