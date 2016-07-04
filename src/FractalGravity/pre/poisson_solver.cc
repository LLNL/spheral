#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void poisson_solver(Group& group, Fractal& fractal,Misc& misc)
  {
    if(group.list_points.size() >= fractal.get_hypre_minimum())
      hypre_solver(group,fractal,misc);
    else
      sor_solver(group,fractal,misc);
  }
}
