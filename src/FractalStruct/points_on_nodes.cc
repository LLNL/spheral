#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void points_on_nodes(Fractal_Memory& mem)
  {
    int lmax=mem.global_level_max;
    vector <int>count;
    mem.p_mess->count_on_node.clear();
    for(int lev=0;lev <= lmax;lev++)
      {
	for(int ni : {0,1})
	  {
	    bool buff=ni > 0;
	    int ok=0;
	    for(auto pg : mem.all_groups[lev])
	      if(pg->get_buffer_group() == buff)
		{
		  ok=1;
		  break;
		}
	    count.push_back(ok);
	    mem.p_mess->count_on_node.push_back(count.back() > 0);
	  }
      }
    vector <int>counts(mem.p_mess->FractalNodes*(lmax+1)*2,0);
    mem.p_mess->my_AllgatherI(count,counts,(lmax+1)*2);
    mem.p_mess->counts_on_nodes.clear();
    mem.p_mess->counts_on_nodes.resize(2*(lmax+1));
    int levelmax=0;
    for(int lev=0;lev <= lmax;lev++)
      {
	for(int ni : {0,1})
	  {
	    for(int FR=0;FR<mem.p_mess->FractalNodes;FR++)
	      {
		int nc=ni+2*lev+2*(lmax+1)*FR;
		if(counts[nc] > 0)
		  levelmax=lev;
		mem.p_mess->counts_on_nodes[ni+2*lev].push_back(counts[nc] > 0);
	      }
	    // mem.p_file->DUMPS << " COUNTP " << mem.steps << " " << lev << " " << ni << " " << mem.p_mess->count_on_node[ni+2*lev] << "\n";
	  }
      }
    // mem.p_file->DUMPS << " levels with good data " << levelmax << " " << mem.global_level_max << "\n";
    //    mem.global_level_max=levelmax;
  }
}
