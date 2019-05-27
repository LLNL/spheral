#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
//
namespace FractalSpace
{
  Misc* Misc::p_misc_instance=NULL;
  bool Misc::get_debug() const
  {
    return debug;
  }
  void Misc::set_debug(bool& d)
  {
    debug=d;
  }
}
