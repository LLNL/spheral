#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
//
namespace FractalSpace
{
  template <class GO_AWAY> void really_clear(vector <GO_AWAY>& die)
  {
    vector <GO_AWAY> YES_REALLY;
    die.swap(YES_REALLY);
  }
}
