#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void poisson_solver(Fractal& fractal,Fractal_Memory& mem,const int& level)
  {
    if(mem.MPIrun)
      {
	hypre_solver(fractal,mem,level);
      }
    else
      {
	for(vector <Group*>::const_iterator group_itr=mem.all_groups[level].begin();
	    group_itr!=mem.all_groups[level].end();group_itr++)
	  {
	    Group& group=**group_itr;
	    sor_solver(group,fractal);
	  }
      }
  }
}
