#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  bool it_is_inside(Point* p_point)
  {
    return p_point != 0 && p_point->get_inside();
  }
}
