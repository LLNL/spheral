#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
//
namespace FractalSpace
{
  template <class GO_AWAY> void really_clear(vector <GO_AWAY>& die)
  {
    die.clear();
    die.shrink_to_fit();
  }
}
namespace FractalSpace
{
  template <class GO_AWAY> void really_resize(vector <GO_AWAY>& die,int howbig)
  {
    die.resize(how_big);
    die.shrink_to_fit();
  }
}
namespace FractalSpace
{
  template <class GO_AWAY> void really_resize2(vector <vector <GO_AWAY> >& die,int howbig)
  {
    die.resize(how_big);
    die.shrink_to_fit();
    for(vector <auto>& d : die)
      d.shrink_to_fit();
  }
}
