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
    vector <bool>done_group;
    int count_sor=0;
    for(vector <Group*>::const_iterator group_itr=mem.all_groups[level].begin();
	group_itr!=mem.all_groups[level].end();group_itr++)
      {
	Group& group=**group_itr;
	FH << "sor 0 " << &group << " " << group.list_points.size() << " " << group.get_buffer_group() << endl;
	if(do_sor || (group.list_points.size() <= mem.min_hypre_group_size && !group.get_buffer_group()))
	  {
	    FH << "sor a " << &group << " " << group.list_points.size() << endl;
	    sor_solver(group,fractal);
	    done_group.push_back(true);
	    FH << "sor b " << &group << " " << group.list_points.size() << " " << count_sor << endl;
	  }
	else
	  done_group.push_back(false);
	count_sor++;
      }
    vector <Point*> p_points_left;
    vector <Point*> p_points_right;
    int counter=0;
    for(vector <Group*>::const_iterator group_itr=mem.all_groups[level].begin();
	group_itr!=mem.all_groups[level].end();group_itr++)
      {
	Group& group=**group_itr;
	FH << "hyp boxes 0 " << &group << " " << group.list_points.size() << " " << counter << " " << done_group[counter] << endl;
	if(!done_group[counter] && !group.get_buffer_group())
	  {
	    FH << "hyp boxes a " << &group << " " << group.list_points.size() << " " << counter << endl;
	    done_group[counter]=true;
	    hypre_little_boxes(group,p_points_left,p_points_right);
	    FH << "hyp boxes b " << &group << " " << group.list_points.size() << " " << counter << endl;
	  }
	counter++;
      }
    FH << " hypre solver self a " << endl;
    if(p_points_left.size() > 0)
      hypre_struct_solver(p_points_left,p_points_right,fractal,mem,level,false);
    FH << " hypre solver self b " << endl;
    p_points_left.clear();
    p_points_right.clear();
    counter=0;
    for(vector <Group*>::const_iterator group_itr=mem.all_groups[level].begin();
	group_itr!=mem.all_groups[level].end();group_itr++)
      {
	Group& group=**group_itr;
	if(!done_group[counter])
	  hypre_little_boxes(group,p_points_left,p_points_right);
	counter++;
      }
    if(p_points_left.size() > 0)
      hypre_struct_solver(p_points_left,p_points_right,fractal,mem,level,true);
  }
}
