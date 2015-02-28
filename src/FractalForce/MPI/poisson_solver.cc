#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void poisson_solver(Fractal& fractal,Fractal_Memory& mem,const int& level)
  {
    if(mem.MPIrun && mem.min_hypre_group_size > 0)
      {
	mem.p_file->FileHypre << " enter MPI Hypre " << level << " " << mem.min_hypre_group_size << endl;
	hypre_solver(fractal,mem,level);
	mem.p_file->FileHypre << "  exit MPI Hypre " << level << " " << mem.min_hypre_group_size << endl;
	return;
      }
    mem.p_file->FileHypre << " Just SOR " << level << endl;
    for(vector <Group*>::const_iterator group_itr=mem.all_groups[level].begin();
	group_itr!=mem.all_groups[level].end();group_itr++)
      {
	Group& group=**group_itr;
	sor_solver(group,fractal);
      }
  }
}
