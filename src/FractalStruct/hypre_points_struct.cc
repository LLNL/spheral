#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void hypre_points_struct(Fractal_Memory& mem,vector <Group*>& groups,
			   vector < vector <Point*> >& hypre_points,bool buffer_groups,int level)
  {
    static int _COUNTER=0;
    mem.hypre_min_node_load=min(mem.hypre_min_node_load,63);
    // ofstream& FHT=mem.p_file->DUMPS;
    vector <int>pos(3);
    vector <int> BOX=mem.BoxesLev[mem.p_mess->FractalRank][level];
    hypre_points.clear();
    for(Group* &pgroup : groups)
      {
 	if(buffer_groups == pgroup->get_buffer_group())
	  {
	    if(!buffer_groups && pgroup->list_points.size() <= mem.hypre_min_node_load)
	      if(mini_solve(mem,pgroup))
		continue;
	    hypre_points.resize(hypre_points.size()+1);
	    for(Point* &p : pgroup->list_points)
	      {
		p->get_pos_point(pos);
		if(p->get_inside() && vector_in_box(pos,BOX))
		  hypre_points.back().push_back(p);
	      }
	    if(!hypre_points.back().empty())
	      sort3_list(hypre_points.back(),0);
	    else
	      hypre_points.resize(hypre_points.size()-1);
	  }
      }
    _COUNTER++;
  }
}
