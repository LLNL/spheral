#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  int findPointByPos(vector <Point*>& plist,Point* psend,ofstream& FF)
  {
    //    FF << "find this " << psend->get_pos_point(0) << " " << psend->get_pos_point(1) << " " << psend->get_pos_point(2) << endl;
    vector <Point*>::iterator found=std::lower_bound(plist.begin(),plist.end(),psend,LesserPoint);
    if(found == plist.end())
      {
	FF << " not found " << endl;
	return -1;
      }
    int label= (found-plist.begin());
    //    FF << "found this " << label << " " << plist[label] << " " << plist[label]->get_pos_point(0) << " " << plist[label]->get_pos_point(1) << " " << plist[label]->get_pos_point(2) << endl;
    return label;
  }
}
