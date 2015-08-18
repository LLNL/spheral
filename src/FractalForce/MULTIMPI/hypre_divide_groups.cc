#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  hypre_divide_groups(vector <Group*>& groups_SOR,vector <Group*>& groups_HYPS,
		      vector <Group*>& groups_HYPM,vector <Group*>& groups,const int& hsize)
  {
    groups_SOR.clear();
    groups_HYPS.clear();
    groups_HYPM.clear();
    if(hsize <= 0)
      {
	groups_SOR=groups;
	return;
      }
    for(vector <Group*>::const_iterator group_itr=groups.begin();group_itr!=groups.end();group_itr++)
      {
	Group* pgroup=*group_itr;
	if(pgroup->list_points.size() <= hsize && !pgroup->get_buffer_group())
	  groups_SOR.push_back(pgroup);
	else if(!pgroup->get_buffer_group())
	  groups_HYPS.push_back(pgroup);
	else
	  groups_HYPM.push_back(pgroup);
      }
  }
}
