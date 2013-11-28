#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void poisson_solver(Fractal& fractal,Fractal_Memory& mem,const int& level)
  {
    bool play_it_safe=false;
    bool do_over=false;
    unsigned int minsize=mem.min_hypre_group_size;
    bool do_sor=minsize <=0 || !mem.MPIrun;
    if(!do_sor)
      {
	mem.p_file->FileHypre << " enter MPI Hypre " << level << " " << mem.min_hypre_group_size << endl;
	if(!play_it_safe)
	  hypre_ij_solver(fractal,mem,level,do_over);
	if(do_over || play_it_safe)
	  hypre_ij_solver_pcg(fractal,mem,level);
	mem.p_file->FileHypre << "  exit MPI Hypre " << level << " " << mem.min_hypre_group_size << endl;
      }
    mem.p_file->FileHypre << " SOR Stuff " << level << endl;
    for(vector <Group*>::const_iterator group_itr=mem.all_groups[level].begin();
	group_itr!=mem.all_groups[level].end();group_itr++)
      {
	Group& group=**group_itr;
	if(do_sor || (group.list_points.size() <= minsize && !group.get_buffer_group()))
	  sor_solver(group,fractal);
      }
  }
}
