#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void points_on_nodes(Fractal_Memory& mem)
  {
    int lmax=mem.global_level_max;
    vector <int>count;
    for(int lev=0;lev <= lmax;lev++)
      {
	bool buff=false;
	for(int ni : {0,1})
	  {
	    int ok=0;
	    for(auto pg : mem.all_groups[lev])
	      if(pg->get_buffer_group() == buff)
		{
		  ok=1;
		  break;
		}
	    count.push_back(ok);
	    buff=true;
	    mem.p_mess->count_on_node.push_back(count.back() > 0);
	  }
      }
    vector <int>counts(mem.p_mess->FractalNodes*(lmax+1)*2);
    mem.p_mess->my_AllgatherI(count,counts,(lmax+1)*2);
    mem.p_mess->counts_on_nodes.clear();
    mem.p_mess->counts_on_nodes.resize(2*(lmax+1));
    int ni=0;
    for(int lev=0;lev <= lmax;lev++)
      for(int FR=0;FR<mem.p_mess->FractalNodes*2;FR++)
	mem.p_mess->counts_on_nodes[lev].push_back(counts[ni++] > 0);
  }
}
