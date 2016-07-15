#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  int fractal_interface(Fractal_Memory* p_fractal_memory,JMO& jmostuff)
  {
    Fractal* p_fractal=0;
    int result=fractal_force_init(p_fractal_memory,p_fractal);
    Fractal_Memory& fractal_memory=*p_fractal_memory;
    Fractal& fractal=*p_fractal;
    populate_fractal(jmostuff,fractal)
    fractal.timing(-2,0);
    fractal.timing(-1,29);
    fractal_force(fractal,fractal_memory);
    fractal.timing(1,29);
    fractal.timing(0,0);
    fractal.timing_lev(0,0);
    populate_JMO(jmostuff,fractal)
    delete p_fractal;
    p_fractal=0;
    return 1;
  }
}
