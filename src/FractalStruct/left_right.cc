#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void left_right(Fractal& frac,vector <double>& pos_left,vector <double>& pos_right)
  {
    pos_left.assign(3,1.0e30);
    pos_right.assign(3,-1.0e30);
    vector <double> pos(3);
    for(int particle=0; particle < frac.get_number_particles(); ++particle)
      {
	frac.particle_list[particle]->get_pos(pos);
	for(int ni=0;ni<3;ni++)
	  {
	    pos_left[ni]=min(pos_left[ni],pos[ni]);
	    pos_right[ni]=max(pos_right[ni],pos[ni]);
	  }
      }
  }	
}
namespace FractalSpace
{
  void left_right(vector <Group*>all_groups,vector <int>& pos_left,vector <int>& pos_right)
  {
    pos_left.assign(3,INT_MAX);
    pos_right.assign(3,INT_MIN);
    vector <int> pos(3);
    // for(vector <Group*>::const_iterator group_itr=all_groups.begin();
    // 	group_itr!=all_groups.end();group_itr++)
    for(auto pg : all_groups)
      {
	// Group& group=**group_itr;
	// for(vector<Point*>::const_iterator point_itr=group.list_points.begin();point_itr !=group.list_points.end();++point_itr)
	for(auto pp : pg->list_points)
	  {
	    // (*point_itr)->get_pos_point(pos);
	    pp->get_pos_point(pos);
	    for(int ni=0;ni<3;ni++)
	      {
		pos_left[ni]=min(pos_left[ni],pos[ni]);
		pos_right[ni]=max(pos_right[ni],pos[ni]);
	      }
	  }
      }
  }	
}
namespace FractalSpace
{
  void left_right(vector <Point*>& all_points,vector <int>& pos_left,vector <int>& pos_right)
  {
    pos_left.assign(3,INT_MAX);
    pos_right.assign(3,INT_MIN);
    vector <int> pos(3);
    // for(vector<Point*>::const_iterator point_itr=all_points.begin();point_itr !=all_points.end();++point_itr)
    for(auto pp : all_points)
      {
	// (*point_itr)->get_pos_point(pos);
	pp->get_pos_point(pos);
	for(int ni=0;ni<3;ni++)
	  {
	    pos_left[ni]=min(pos_left[ni],pos[ni]);
	    pos_right[ni]=max(pos_right[ni],pos[ni]);
	  }
      }
  }	
}
namespace FractalSpace
{
  void left_right(vector <Point*>& all_points,vector <int>& pos_left,vector <int>& pos_right,const int& wrap)
  {
    pos_left.assign(3,INT_MAX);
    pos_right.assign(3,INT_MIN);
    vector <int> pos(3);
    // for(vector<Point*>::const_iterator point_itr=all_points.begin();point_itr !=all_points.end();++point_itr)
    for(auto pp : all_points)
      {
	// (*point_itr)->get_pos_point(pos);
	pp->get_pos_point(pos);
	for(int ni=0;ni<3;ni++)
	  {
	    pos[ni]=(pos[ni]+wrap) % wrap;
	    pos_left[ni]=min(pos_left[ni],pos[ni]);
	    pos_right[ni]=max(pos_right[ni],pos[ni]);
	  }
      }
  }
}
namespace FractalSpace
{
  void left_right(deque <Point*>& all_points,vector <int>& pos_left,vector <int>& pos_right)
  {
    pos_left.assign(3,INT_MAX);
    pos_right.assign(3,INT_MIN);
    vector <int> pos(3);
    // for(vector<Point*>::const_iterator point_itr=all_points.begin();point_itr !=all_points.end();++point_itr)
    for(auto pp : all_points)
      {
	// (*point_itr)->get_pos_point(pos);
	pp->get_pos_point(pos);
	for(int ni=0;ni<3;ni++)
	  {
	    pos_left[ni]=min(pos_left[ni],pos[ni]);
	    pos_right[ni]=max(pos_right[ni],pos[ni]);
	  }
      }
  }	
}
namespace FractalSpace
{
  void left_right(deque <Point*>& all_points,vector <int>& pos_left,vector <int>& pos_right,const int& wrap)
  {
    pos_left.assign(3,INT_MAX);
    pos_right.assign(3,INT_MIN);
    vector <int> pos(3);
    // for(vector<Point*>::const_iterator point_itr=all_points.begin();point_itr !=all_points.end();++point_itr)
    for(auto pp : all_points)
      {
	// (*point_itr)->get_pos_point(pos);
	pp->get_pos_point(pos);
	for(int ni=0;ni<3;ni++)
	  {
	    pos[ni]=(pos[ni]+wrap) % wrap;
	    pos_left[ni]=min(pos_left[ni],pos[ni]);
	    pos_right[ni]=max(pos_right[ni],pos[ni]);
	  }
      }
  }	
}
