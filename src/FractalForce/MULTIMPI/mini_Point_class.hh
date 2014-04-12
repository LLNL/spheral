#ifndef _mini_Point_Defined_
#define _mini_Point_Defined_
namespace FractalSpace
{
  class mini_Point{
  public:
    bool realpoint;
    bool it_is_a_point;
    Point* pmyself;
    mini_Point():
      realpoint(false),
      it_is_a_point(false),
      pmyself(0)
    {}
    ~mini_Point()
    {}
  };
}
#endif
