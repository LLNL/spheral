#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
//
namespace FractalSpace
{
  void tree_start(Group& group, Fractal& fractal,Fractal_Memory& mem,Misc& misc)
  {
    bool halo=fractal.get_halo();
    double a_grid_length=(double)fractal.get_grid_length();
    vector <int>Boxu(6);
    vector <Point*>ud(6);
    Point* haha=0;
    double t0,t1,t2,t3;
    t0=clock();
    cout << "enter treestart" << endl;
    int point_counter=0;
    mem.total_points_generated=0;
    mem.total_points_used=0;
    Point* p_point=0;
    Point* new_points=0;
    int new_counter=0;
    cout << " generating points in treestart" << endl;
    vector <int> Box;
    fractal.getBox(Box);
    vector <int> BoxReal;
    fractal.getBoxReal(BoxReal);
    int volume=(Box[1]-Box[0]+1);
    volume*=(Box[3]-Box[2]+1);
    volume*=(Box[5]-Box[4]+1);
    new_points=new (nothrow) Point[volume];
    assert(new_points);
    mem.total_points_generated=volume;
    group.list_new_points.push_back(new_points);
    group.list_points.reserve(volume);
    cout << " generated points in treestart" << endl;
    cout << " Box ";
    cout << Box[0] << " ";
    cout << Box[1] << " ";
    cout << Box[2] << " ";
    cout << Box[3] << " ";
    cout << Box[4] << " ";
    cout << Box[5] << endl;
    cout << " BoxReal ";
    cout << BoxReal[0] << " ";
    cout << BoxReal[1] << " ";
    cout << BoxReal[2] << " ";
    cout << BoxReal[3] << " ";
    cout << BoxReal[4] << " ";
    cout << BoxReal[5] << endl;
    //
    for(int grid_x=Box[0];grid_x <= Box[1];grid_x++)
      {
	bool edge_x= grid_x == Box[0] || grid_x == Box[1];
	bool buff_x= grid_x < BoxReal[0] || grid_x > BoxReal[1];
	for(int grid_y=Box[2];grid_y <= Box[3];grid_y++)
	  {
	    bool edge_y= grid_y == Box[2] || grid_y == Box[3];
	    bool buff_y= grid_y < BoxReal[2] || grid_y > BoxReal[3];
	    for(int grid_z=Box[4];grid_z <= Box[5];grid_z++)
	      {
		bool edge_z= grid_z == Box[4] || grid_z == Box[5];
		bool buff_z= grid_z < BoxReal[4] || grid_z > BoxReal[5];
		bool edge=edge_x || edge_y || edge_z;
		bool buff=buff_x || buff_y || buff_z;
		p_point=&new_points[new_counter];
		new_counter++;
		mem.total_points_used++;
		Point& point=*p_point;
		group.list_points.push_back(p_point);
		point.set_point_to_number(point_counter);
		point_counter++;
		point.set_point_pointer(p_point);
		point.set_real_pointer(0);
		point.set_pos_point(grid_x*misc.zoom,grid_y*misc.zoom,grid_z*misc.zoom);
		point.set_inside(!edge);
		if(edge) point.set_real_pointer(-1);
		point.set_p_in_group(misc.p_group_0);
		point.set_p_in_high_group(0);
		point.set_buffer_surface_point(buff,edge);
	      }
	  }
      }
    cout << "aa " << point_counter << endl;
    t1=clock();
    cout << "treestart enter second part" << endl;
    int grid=0;
    for(int grid_x=Box[0];grid_x <= Box[1];grid_x++)
      {
	for(int grid_y=Box[2];grid_y <= Box[3];grid_y++)
	  {
	    for(int grid_z=Box[4];grid_z <= Box[5];grid_z++)
	      {
		ud.assign(6,haha);
		where_6(grid_x,grid_y,grid_z,Box,Boxu);
		for(int ni=0;ni<6;ni++)
		  {
		    if(Boxu[ni] >= 0)
		      ud[ni]=group.list_points[Boxu[ni]];
		    else
		      ud[ni]=0;
		  }
		Point* p_point=group.list_points[grid];
		p_point->set_point_ud(ud);
		grid++;
	      }
	  }
      }
    cout << "bb " << grid << endl;
    cout << "ccc " << endl;
    if(misc.get_debug())
      {
	check_for_edge_trouble(fractal);
	cout << "checked for edge" << endl;
      }
    t2=clock();
    //
    int moat=max(max(1,fractal.get_moat_0()),fractal.get_padding());
    vector <int>Boxm;
    Boxm=Box;
    Boxm[0]=Box[0]+moat;
    Boxm[1]=Box[1]-moat;
    Boxm[2]=Box[2]+moat;
    Boxm[3]=Box[0]-moat;
    Boxm[4]=Box[4]+moat;
    Boxm[5]=Box[5]-moat;
    double radmin2=0.0;
    if(halo) radmin2=pow(0.5-(double)moat/(double)fractal.get_grid_length(),2);
    vector <double> pos(3);
    for(int particle=0; particle < fractal.get_number_particles(); ++particle)
      {
	Particle* p=fractal.particle_list[particle];
	p->set_p_highest_level_group(misc.p_group_0);
	p->set_highest_level(0);
	p->get_pos(pos);
	int grid_x=(int)floor(pos[0]*a_grid_length);
	int grid_y=(int)floor(pos[1]*a_grid_length);
	int grid_z=(int)floor(pos[2]*a_grid_length);
	bool do_it=
	  grid_x >= Box[0] && grid_x < Box[1] &&
	  grid_y >= Box[2] && grid_y < Box[3] &&
	  grid_z >= Box[4] && grid_z < Box[5];
	bool do_it_ins=
	  grid_x >= BoxReal[0] && grid_x <= BoxReal[1] &&
	  grid_y >= BoxReal[2] && grid_y <= BoxReal[3] &&
	  grid_z >= BoxReal[4] && grid_z <= BoxReal[5];
	if(halo && do_it)
	  do_it=p->get_r2(0.5,0.5,0.5) < radmin2;
	if(halo && do_it_ins)
	  do_it_ins=p->get_r2(0.5,0.5,0.5) < radmin2;
	p->set_buffer_particle(!do_it_ins);
	if(do_it)
	  {
	    int grid=where_1(grid_x,grid_y,grid_z,Box);
	    //	    cout << "grid " << grid_x << " "  << grid_y << " "  << grid_z << " " << grid << endl;
	    //	    cout << pos[0] << " "  << pos[1] << " "  << pos[2] << endl;
	    assert(grid >= 0);
	    Point* p_point=group.list_points[grid];
	    p_point->list_particles.push_back(p);
	  }
	else
	  p->set_p_highest_level_group(0);
      }
    cout << "made it here " << endl;
    cout << "end tree " << " " << group.list_points.size() << endl;
    if(fractal.get_level_max() > 0)
      group.set_number_high_groups(-1);
    else
      group.set_number_high_groups(0);
    cout << "exit treestart" << endl;
    t3=clock();
    cout << "tree time " << t1-t0 << " " << t2-t1 << " " << t3-t2 << endl;
    /***
    vector <int> posp(3);
    for(vector<Point*>::const_iterator point_itr=group.list_points.begin();point_itr !=group.list_points.end();++point_itr)
      {
	Point& point=**point_itr;
	if(point.list_particles.empty()) continue;
	point.get_pos_point(posp);
	cout << " test point " << posp[0] << " " << posp[1] << " " << posp[2] << endl;
	for(vector<Particle*>::const_iterator particle_itr=point.list_particles.begin();particle_itr !=point.list_particles.end();++particle_itr)
	  {
	    Particle& particle=**particle_itr;
	    particle.get_pos(pos);
	    cout << " test particle " << pos[0] << " " << pos[1] << " " << pos[2] << endl;
	  }
      }
    */
  }
  //
  void where_6(const int& i,const int& j,const int& k,vector <int>& Box,vector <int>& Boxu)
  {
    Boxu.assign(6,-1);
    Boxu[0]=where_1(i-1,j,k,Box);
    Boxu[1]=where_1(i+1,j,k,Box);
    Boxu[2]=where_1(i,j-1,k,Box);
    Boxu[3]=where_1(i,j+1,k,Box);
    Boxu[4]=where_1(i,j,k-1,Box);
    Boxu[5]=where_1(i,j,k+1,Box);
  }
  int where_1(const int& i,const int& j,const int& k,vector <int>& Box)
  {
    if(i >= Box[0] && i <= Box[1] && j >= Box[2] && j <= Box[3] && k >= Box[4] && k <= Box[5])
      return (k-Box[4])+((j-Box[2])+(i-Box[0])*(Box[3]-Box[2]+1))*(Box[5]-Box[4]+1);
    return -1;
  }
}
