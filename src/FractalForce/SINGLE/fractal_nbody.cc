#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
int main()
{
  using namespace FractalSpace;
  cout << "starting out " << endl;
  Fractal_Memory* p_fractal_memory= 0;
  Fractal* p_fractal=0;
  int result=fractal_gravity_wrapper(p_fractal_memory,p_fractal);
  delete p_fractal_memory;
  p_fractal_memory=0;
  delete p_fractal;
  p_fractal=0;
  return result;
}
