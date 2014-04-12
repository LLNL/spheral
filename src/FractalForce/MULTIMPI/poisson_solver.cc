#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void poisson_solver(Fractal& fractal,Fractal_Memory& mem,const int& level)
  {
    FILE* PFH=mem.p_file->PFHypre;
    unsigned int minsize=mem.min_hypre_group_size;
    fprintf(PFH," enter MPI Hypre %d %d \n",level, mem.min_hypre_group_size);
    hypre_ij_solver(fractal,mem,level);
    fprintf(PFH," exit MPI Hypre %d %d \n",level, mem.min_hypre_group_size);
    fprintf(PFH," SOR Stuff %d \n",level);
    for(vector <Group*>::const_iterator group_itr=mem.all_groups[level].begin();
	group_itr!=mem.all_groups[level].end();group_itr++)
      {
	Group& group=**group_itr;
	if(group.list_points.size() <= minsize && !group.get_buffer_group())
	  sor_solver(group,fractal);
      }
  }
}
