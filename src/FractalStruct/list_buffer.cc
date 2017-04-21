#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  //
  void list_buffer(Point& point,const int& corner)
  {
    assert(corner >=0 && corner < 8);
    Point* black_knight=&point;
    // for(int ni=0;ni < 7;++ni)
    //   {
    // 	black_knight=black_knight->get_point_ud_0(Point::order[corner][ni],1);
    // 	black_knight->set_it_is_high(true);
    //   }
    for (auto pN : Point::order[corner])
      {
    	black_knight=black_knight->get_point_ud_0(pN,1);
    	black_knight->set_it_is_high(true);
      }
  }
}
