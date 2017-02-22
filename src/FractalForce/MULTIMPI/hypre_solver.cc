#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void hypre_solver(Fractal& fractal,Fractal_Memory& mem,const int& level)
  {
    ofstream& FH=mem.p_file->FileHypre;
    bool do_sor=mem.min_hypre_group_size <= 0;
    FH << "hypre solver a " << level << " " << do_sor << " " << mem.min_hypre_group_size << endl;
    vector <Group*>groups_SOR;
    vector <Group*>groups_HYPS;
    vector <Group*>groups_HYPM;
    if(do_sor)
      groups_SOR=mem.all_groups[level];
    else
      hypre_divide_groups(groups_SOR,groups_HYPS,groups_HYPM,mem.all_groups[level],mem.hypre_group_size);
    end if
      FH << " HYP sizes " << groups_SOR.size() << " " << groups_HYPS.size() << " " << groups_HYPM.size() << endl;
    for(vector <Group*>::const_iterator group_itr=groups_SOR.begin();
	group_itr!=mem.groups_SOR.end();group_itr++)
      {
	Group& group=**group_itr;
	FH << " sor a " << &group << " " << group.list_points.size() << " " << group.get_buffer_group() << endl;
	sor_solver(group,fractal);
	FH << " sor b " << &group << " " << group.list_points.size() << " " << group.get_buffer_group() << endl;
      }
    vector <Point*> p_points_left;
    vector <Point*> p_points_right;
    if(groups_HYPS.size() > 0)
      {
	for(vector <Group*>::const_iterator group_itr=groups_HYPS.begin();
	    group_itr!=mem.groups_HYPS.end();group_itr++)
	  {
	    Group& group=**group_itr;
	    FH << " HYPS a " << &group << " " << group.list_points.size() << " " << group.get_buffer_group() << endl;
	    hypre_big_boxes(fractal,group,p_points_left,p_points_right);
	    FH << " HYPS b " << &group << " " << group.list_points.size() << " " << group.get_buffer_group() << endl;
	  }
	FH << " HYPS C " << p_points_left.size() << endl;
	hypre_struct_solver(p_points_left,p_points_right,fractal,mem,level,false);
	FH << " HYPS D " << p_points_left.size() << endl;
      }
    p_points_left.clear();
    p_points_right.clear();
    if(groups_HYPM.size() > 0)
      {
	for(vector <Group*>::const_iterator group_itr=groups_HYPM.begin();
	    group_itr!=mem.groups_HYPM.end();group_itr++)
	  {
	    Group& group=**group_itr;
	    FH << " HYPM a " << &group << " " << group.list_points.size() << " " << group.get_buffer_group() << endl;
	    hypre_big_boxes(fractal,group,p_points_left,p_points_right);
	    FH << " HYPM B " << &group << " " << group.list_points.size() << " " << group.get_buffer_group() << endl;
	  }
	FH << " HYPM C " << p_points_left.size() << endl;
	hypre_struct_solver(p_points_left,p_points_right,fractal,mem,level,true);
	FH << " HYPM C " << p_points_left.size() << endl;
      }
  }
}
