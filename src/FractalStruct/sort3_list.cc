#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  int Misc::dim0=0;
  int Misc::dim1=1;
  int Misc::dim2=2;
  //
  void sort3_list(Group& group,int what)
  {
    Misc::dim2=(what+2) % 3;
    Misc::dim1=(what+1) % 3;
    Misc::dim0=(what+0) % 3;
    sort(group.list_points.begin(),group.list_points.end(),LesserPoint);
  }
  void sort3_list(vector <Point*>& list_points,int what)
  {
    Misc::dim2=(what+2) % 3;
    Misc::dim1=(what+1) % 3;
    Misc::dim0=(what+0) % 3;
    sort(list_points.begin(),list_points.end(),LesserPoint);
  }
  void sort3_list(deque <Point*>& list_points,int what)
  {
    Misc::dim2=(what+2) % 3;
    Misc::dim1=(what+1) % 3;
    Misc::dim0=(what+0) % 3;
    sort(list_points.begin(),list_points.end(),LesserPoint);
  }
  void sort3_list(list <Point*>& list_points,int what)
  {
    Misc::dim2=(what+2) % 3;
    Misc::dim1=(what+1) % 3;
    Misc::dim0=(what+0) % 3;
    list_points.sort(LesserPoint);
  }
  bool LesserPoint(Point* p1,Point* p2)
  {
    int dz(p1->get_pos_point(Misc::dim2)-p2->get_pos_point(Misc::dim2));
    if(dz != 0)
      return dz < 0;
    int dy(p1->get_pos_point(Misc::dim1)-p2->get_pos_point(Misc::dim1));
    if(dy != 0)
      return dy < 0;
    return (p1->get_pos_point(Misc::dim0)-p2->get_pos_point(Misc::dim0) < 0);
  }
  bool LesserPointA(Point* p1,Point* p2)
  {
    //    FF << "Lesser A " << p1 << " " << p1->get_pos_point(0) << " " << p1->get_pos_point(1) << " " << p1->get_pos_point(2) << "\n";
    //    FF << "Lesser B " << p2 << " " << p2->get_pos_point(0) << " " << p2->get_pos_point(1) << " " << p2->get_pos_point(2) << "\n";
    int dz=p1->get_pos_point(Misc::dim2)-p2->get_pos_point(Misc::dim2);
    if(dz != 0)
      return dz < 0;
    int dy=p1->get_pos_point(Misc::dim1)-p2->get_pos_point(Misc::dim1);
    if(dy != 0)
      return dy < 0;
    int dx=p1->get_pos_point(Misc::dim0)-p2->get_pos_point(Misc::dim0);
    return dx < 0;
  }
}
