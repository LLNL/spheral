#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void fractal_force_init(Fractal_Memory* PFM)
  {
    PFM->calc_FractalNodes();
    PFM->calc_Buffers_and_more();
    PFM->calc_RealBoxes();
  }
}
 
