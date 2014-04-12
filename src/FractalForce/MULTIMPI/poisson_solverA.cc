#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void poisson_solver(Fractal& fractal,Fractal_Memory& mem,const int& level)
  {
    // bool buffer_only (false everything, true buffer groups only)
    FILE* PFH=mem.p_file->PFHypre;
    int m_size=mem.min_hypre_group_size;
    fprintf(PFH," enter MPI Hypre a %d %d \n",level, m_size);
    //    hypre_ij_solver(fractal,mem,level,true);
    //    hypre_ij_solver(fractal,mem,level,false);
    fprintf(PFH," enter MPI Hypre b %d %d \n",level, m_size);
    hypre_ij_solver_selfie(fractal,mem,level);
    fprintf(PFH," exit MPI Hypre b %d %d \n",level, m_size);
    mem.p_mess->Full_Stop();
    fprintf(PFH," SOR Stuff enter %d \n",level);
    int ngroups=0;
    for(vector <Group*>::const_iterator group_itr=mem.all_groups[level].begin();
	group_itr!=mem.all_groups[level].end();group_itr++)
      {
	Group& group=**group_itr;
	if(group.list_points.size() <= m_size && !group.get_buffer_group())
	  {
	    sor_solver(group,fractal);
	    ngroups++;
	  }
      }
    fprintf(PFH," SOR Stuff exit %d %d \n",level,ngroups);
  }
}
