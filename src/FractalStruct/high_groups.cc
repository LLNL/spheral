#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void high_groups(Group& group)
  {
    int n_high_p=group.get_number_high_points();
    vector <Group*> family(n_high_p);
    Group* p_group=&group;
    for(int n=0; n < n_high_p;++n)
      {
	if(group.head_number[n]==n)
	  {
	    Group* p_high_group=new Group;
	    assert(p_high_group);
	    Group& high_group=*p_high_group;
	    group.list_high_groups.push_back(p_high_group);
	    high_group.set_p_generated_from_group(p_group);
	    high_group.list_high_points.push_back(group.list_high[n]);
	    family[n]=p_high_group;
	    group.list_high[n]->set_p_in_high_group(family[n]);
	  }
      }
    for(int n=0; n < n_high_p;++n)
      {
	if(group.head_number[n]!=n)
	  {
	    Group* p_high_g=family[group.head_number[n]];
	    p_high_g->list_high_points.push_back(group.list_high[n]);
	    group.list_high[n]->set_p_in_high_group(p_high_g);
	  }
      }
    clean_vector(group.head_number);
    clean_vector(group.list_high);
  }
}
