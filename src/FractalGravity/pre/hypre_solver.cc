#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void hypre_solver(Group& group, Fractal& fractal,Misc& misc)
  {
    vector < vector <int> > boxab;
    boxes(fractal,group,boxab);
  }
}
