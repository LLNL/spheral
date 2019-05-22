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
  template <class T> bool vector_in_box(vector <T>& xvec,vector <T>& box)
  {
    return (xvec[0] <= box[1] && xvec[0] >= box[0] && 
	    xvec[1] <= box[3] && xvec[1] >= box[2] &&
	    xvec[2] <= box[5] && xvec[2] >= box[4]);
  }
}
namespace FractalSpace
{
  template bool vector_in_box(vector <int>& xvec,vector <int>& box);
  template bool vector_in_box(vector <double>& xvec,vector <double>& box);
}
namespace FractalSpace
{
  template <class T> bool vector_in_box(vector <T>& xvec,vector <T>& box,T width,bool periodic)
  {
    const int many= periodic ? 1:0;
    for(int mz=-many;mz <= many;mz++)
      {
	vector <T>XVECA(xvec);
	XVECA[2]=xvec[2]+mz*width;
	for(int my=-many;my <= many;my++)
	  {
	    XVECA[1]=xvec[1]+my*width;
	    for(int mx=-many;mx <= many;mx++)
	      {
		XVECA[0]=xvec[0]+mx*width;
		if(vector_in_box(XVECA,box))
		  {
		    xvec=XVECA;
		    return true;
		  }
	      }
	  }
      }
    return false;
  }
}
namespace FractalSpace
{
  template bool vector_in_box(vector <int>& xvec,vector <int>& box,int width,bool periodic);
  template bool vector_in_box(vector <double>& xvec,vector <double>& box,double width,bool periodic);
}
namespace FractalSpace
{
  template <class T> bool vector_in_box(array <T,3>& xvec,vector <T>& box)
  {
    return (xvec[0] <= box[1] && xvec[0] >= box[0] && 
	    xvec[1] <= box[3] && xvec[1] >= box[2] &&
	    xvec[2] <= box[5] && xvec[2] >= box[4]);
  }
}
namespace FractalSpace
{
  template bool vector_in_box(array <int,3>& xvec,vector <int>& box);
  template bool vector_in_box(array <double,3>& xvec,vector <double>& box);
}
namespace FractalSpace
{
  bool vector_in_box(Point* p,vector <int>& box)
  {
    vector <int>xvec(3);
    p->get_pos_point(xvec);
    return (xvec[0] <= box[1] && xvec[0] >= box[0] && 
	    xvec[1] <= box[3] && xvec[1] >= box[2] &&
	    xvec[2] <= box[5] && xvec[2] >= box[4]);
  }
}
namespace FractalSpace
{
  bool vector_in_box(const Point& p,vector <int>& box)
  {
    vector <int>xvec(3);
    p.get_pos_point(xvec);
    return (xvec[0] <= box[1] && xvec[0] >= box[0] && 
	    xvec[1] <= box[3] && xvec[1] >= box[2] &&
	    xvec[2] <= box[5] && xvec[2] >= box[4]);
  }
}
// namespace FractalSpace
// {
//   template <class T> T number_in_box(vector <T>& xvec,vector <T>& box)
//   {
//     return (xvec[0] <= box[1] && xvec[0] >= box[0] && 
// 	    xvec[1] <= box[3] && xvec[1] >= box[2] &&
// 	    xvec[2] <= box[5] && xvec[2] >= box[4]);
//   }
// }
// namespace FractalSpace
// {
//   template int number_in_box(vector <int>& xvec,vector <int>& box);
//   template double number_in_box(vector <double>& xvec,vector <double>& box);
// }
namespace FractalSpace
{
  template <class T> bool overlap_boxes(vector <T>& boxa,vector <T>& boxb)
  {
    return (boxa[0] <= boxb[1] && boxa[1] >= boxb[0] && 
	    boxa[2] <= boxb[3] && boxa[3] >= boxb[2] &&
	    boxa[4] <= boxb[5] && boxa[5] >= boxb[4]);
  }
}
namespace FractalSpace
{
  template bool overlap_boxes(vector <int>& boxa,vector <int>& boxb);
  template bool overlap_boxes(vector <double>& boxa,vector <double>& boxb);
}
namespace FractalSpace
{
  template <class T> bool overlap_boxes(vector <T>& boxa,vector <T>& boxb,T width,bool periodic)
  {
    const int many= periodic ? 1:0;
    //    cerr << " TWD " << many << " " << periodic << "\n";
    // *Point::p_FILE << " TWD " << many << " " << periodic << "\n";
    for(int mz=-many;mz <= many;mz++)
      {
	vector <T>BOXA(boxa);
	BOXA[4]=boxa[4]+mz*width;
	BOXA[5]=boxa[5]+mz*width;
	for(int my=-many;my <= many;my++)
	  {
	    BOXA[2]=boxa[2]+my*width;
	    BOXA[3]=boxa[3]+my*width;
	    for(int mx=-many;mx <= many;mx++)
	      {
		BOXA[0]=boxa[0]+mx*width;
		BOXA[1]=boxa[1]+mx*width;
		// *Point::p_FILE << " TWC " << mx << " " << my << " " << mz << " ";
		// *Point::p_FILE << BOXA[0] << " " << BOXA[1] << " " << BOXA[2] << " " << BOXA[3] << " " << BOXA[4] << " " << BOXA[5] << " ";
		// *Point::p_FILE << boxb[0] << " " << boxb[1] << " " << boxb[2] << " " << boxb[3] << " " << boxb[4] << " " << boxb[5] << "\n";
		if(overlap_boxes(BOXA,boxb))
		  return true;
	      }
	  }
      }
    return false;
  }
}
namespace FractalSpace
{
  template bool overlap_boxes(vector <int>& boxa,vector <int>& boxb,int width,bool periodic);
  template bool overlap_boxes(vector <double>& boxa,vector <double>& boxb,double width,bool periodic);
}
///
namespace FractalSpace
{
  bool group_in_box(Group* pgroup,vector<int>& BOX)
  {
    vector<int>miny;
    vector<int>maxy;
    pgroup->get_miny_maxy(miny,maxy);
    return overlap(miny,maxy,BOX);
  }
}
