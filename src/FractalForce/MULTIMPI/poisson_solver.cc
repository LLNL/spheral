#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void poisson_solver(Fractal& fractal,Fractal_Memory& mem,const int& level)
  {
    double gcells=mem.grid_length;
    gcells=pow(gcells,3);
    double FNO=mem.FractalNodes;
    mem.hypre_max_node_load=4.0*gcells/FNO;
    // bool buffer_only (false everything, true buffer groups only)
    bool buffer_only=fractal.get_periodic();
    //    buffer_only=true;
    buffer_only=false;
    FILE* PFH=mem.p_file->PFHypre;
    ofstream& FHT=mem.p_file->DUMPS;
    int m_size=mem.min_hypre_group_size;
    fprintf(PFH," enter MPI Hypre a %d %d \n",level, m_size);
    hypre_ij_solver(fractal,mem,level,buffer_only);
    //    hypre_ij_solver(fractal,mem,level,false);
    if(buffer_only)
      {
	fprintf(PFH," enter MPI Hypre selfie %d %d \n",level, m_size);
	hypre_ij_solver_selfie(fractal,mem,level);
	fprintf(PFH," exit MPI Hypre selfie %d %d \n",level, m_size);
      }
    //    mem.p_mess->Full_Stop();
    double SOR_time=-mem.p_mess->Clock();
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
    SOR_time+=mem.p_mess->Clock();
    FHT << " S" << mem.steps << "S " << "L" << level << "L" << "\t" << SOR_time << "\t" << "SOR   Time       " << "\n";
    fprintf(PFH," SOR Stuff exit %d %d \n",level,ngroups);
  }
}
