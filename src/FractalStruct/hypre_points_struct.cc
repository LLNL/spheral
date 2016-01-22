#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void hypre_points_struct(Fractal_Memory& mem,vector <Group*>& groups,
			   vector < vector <Point*> >& hypre_points,bool buffer_groups,int level)
  {
    vector <int>pos(3);
    vector <int> BOX=mem.BoxesLev[mem.p_mess->FractalRank][level];
    for(vector <Group*>::const_iterator group_itr=groups.begin();group_itr!=groups.end();group_itr++)
      {
	Group* pgroup=*group_itr;
	if(buffer_groups == pgroup->get_buffer_group())
	  {
	    hypre_points.resize(hypre_points.size()+1);
	    for(vector<Point*>::const_iterator point_itr=pgroup->list_points.begin();point_itr !=pgroup->list_points.end();++point_itr)
	      {
		Point* p=*point_itr;
		p->get_pos_point(pos);
		if(p->get_inside() && vector_in_box(pos,BOX))
		  (hypre_points.back()).push_back(p);
	      }
	    if(!hypre_points.back().empty())
	      sort3_list(hypre_points.back(),0);
	  }
      }
  }
}
