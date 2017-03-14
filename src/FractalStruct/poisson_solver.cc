#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void poisson_solver(Fractal& fractal,Fractal_Memory& mem,const int& level)
  {
//     mem.hypre_max_node_load=50000;
//     ofstream& FHT=mem.p_file->DUMPS;
    //    int m_size=mem.min_hypre_group_size;
    hypre_ij_solver(fractal,mem,level,true);
    hypre_ij_solver(fractal,mem,level,false);
//     double SOR_time=-mem.p_mess->Clock();
//     int ngroups=0;
//     for(vector <Group*>::const_iterator group_itr=mem.all_groups[level].begin();
// 	group_itr!=mem.all_groups[level].end();group_itr++)
//       {
// 	Group& group=**group_itr;
// 	if(group.list_points.size() <= mem.min_hypre_group_size && !group.get_buffer_group())
// 	  {
// 	    sor_solver(group,fractal);
// 	    ngroups++;
// 	  }
//       }
//     SOR_time+=mem.p_mess->Clock();
//     FHT << " S" << mem.steps << "S " << "L" << level << "L" << "\t" << SOR_time << "\t" << "SOR   Time       " << "\n";
  }
}
