#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
//
namespace FractalSpace
{
  void tree_start(Group& group, Fractal& fractal,Fractal_Memory& mem,Misc& misc)
  {
    double a_grid_length=(double)fractal.get_grid_length();
    int lengthm1=fractal.get_grid_length()-1;
    vector <int>Boxu(6);
    vector <Point*>ud(6);
    double t0,t1,t2,t3;
    t0=clock();
    cout << "enter treestart" << endl;
    bool period=fractal.get_periodic();
    vector <bool> Periods(3);
    fractal.get_Periods(Periods);
    int point_counter=0;
    mem.total_points_generated=0;
    mem.total_points_used=0;
    Point* p_point=0;
    Point* new_points=0;
    cout << " generating points in treestart" << endl;
    vector <int> Box;
    fractal.getBox(Box);
    cout << "Box ";
    cout << Box[0] << " ";
    cout << Box[1] << " ";
    cout << Box[2] << " ";
    cout << Box[3] << " ";
    cout << Box[4] << " ";
    cout << Box[5] << " ";
    cout << Periods[0] << Periods[1] << Periods[2];
    cout << endl;
    bool MPIrun=mem.MPIrun;
    group.set_buffer_group(MPIrun);
    if(!MPIrun)
      {
	assert(Box[0] == 0);
	assert(Box[2] == 0);
	assert(Box[4] == 0);
	assert(Box[1] == lengthm1);
	assert(Box[3] == lengthm1);
	assert(Box[5] == lengthm1);
      }
    else
      {
	assert(Box[0] >= 0);
	assert(Box[2] >= 0);
	assert(Box[4] >= 0);
	assert(Box[1] <= lengthm1);
	assert(Box[3] <= lengthm1);
	assert(Box[5] <= lengthm1);
      }
    vector <int> Buffer(6);
    fractal.getBuffer(Buffer);
    vector <int> BBox(6);
    fractal.getBBox(BBox);
    vector <int> PBox(6);
    fractal.getPBox(PBox);
    vector <int> PBoxLength(3);
    fractal.getPBoxLength(PBoxLength);
    int volume=PBoxLength[0]*PBoxLength[1]*PBoxLength[2];
    new_points=new (nothrow) Point[volume];
    assert(new_points);
    mem.total_points_generated=volume;
    group.set_group_number(0);
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
    cout << " Buffer ";
    cout << Buffer[0] << " ";
    cout << Buffer[1] << " ";
    cout << Buffer[2] << " ";
    cout << Buffer[3] << " ";
    cout << Buffer[4] << " ";
    cout << Buffer[5] << endl;
    //
    vector <int>grid(3);
    bool inside;
    bool edge;
    bool buff;
    bool pass;
    int new_counter=0;
    for(int gridx=PBox[0];gridx <= PBox[1];gridx++)
      {
	grid[0]=gridx;
	for(int gridy=PBox[2];gridy <= PBox[3];gridy++)
	  {
	    grid[1]=gridy;
	    for(int gridz=PBox[4];gridz <= PBox[5];gridz++)
	      {
		grid[2]=gridz;
		p_point=&new_points[new_counter];
		new_counter++;
		mem.total_points_used++;
		Point& point=*p_point;
		group.list_points.push_back(p_point);
		point.set_point_to_number(point_counter);
		point.set_point_pointer(p_point);
		point.set_real_pointer(0);
		point.set_pos_point(gridx*misc.zoom,gridy*misc.zoom,gridz*misc.zoom);
		point.set_p_in_group(misc.p_group_0);
		point.set_p_in_high_group(0);
		fractal.inside_edge_buffer_pass(grid,inside,edge,buff,pass);
		point.set_inside(inside);
		point.set_edge_point(edge);
		point.set_buffer_point(buff);
		point.set_buffer_point(pass);
		point_counter++;
	      }
	  }
      }
    cout << "aa " << point_counter << endl;
    t1=clock();
    cout << "treestart enter second part" << endl;
    int gridxyz=0;
    for(int grid_x=PBox[0];grid_x <= PBox[1];grid_x++)
      {
	for(int grid_y=PBox[2];grid_y <= PBox[3];grid_y++)
	  {
	    for(int grid_z=PBox[4];grid_z <= PBox[5];grid_z++)
	      {
		ud.assign(6,Point::nothing);
		fractal.where_6(grid_x,grid_y,grid_z,Boxu);
		for(int ni=0;ni<6;ni++)
		  {
		    if(Boxu[ni] >= 0)
		      ud[ni]=group.list_points[Boxu[ni]];
		    else
		      ud[ni]=0;
		  }
		Point* p_point=group.list_points[gridxyz];
		p_point->set_point_ud(ud);
		gridxyz++;
	      }
	  }
      }
    cout << "bb " << gridxyz << endl;
    if(misc.get_debug())
      {
	check_for_edge_trouble(fractal);
	cout << "checked for edge" << endl;
      }
    t2=clock();
    //
    int moat=max(max(1,fractal.get_moat_0()),fractal.get_padding());
    int glm=fractal.get_grid_length()-moat;
    vector <double> pos(3);
    for(int particle=0; particle < fractal.get_number_particles(); ++particle)
      {
	Particle* p=fractal.particle_list[particle];
	p->set_p_highest_level_group(misc.p_group_0);
	p->set_highest_level(0);
	p->get_pos(pos);
	int grid_x=static_cast<int>(floor(pos[0]*a_grid_length));
	int grid_y=static_cast<int>(floor(pos[1]*a_grid_length));
	int grid_z=static_cast<int>(floor(pos[2]*a_grid_length));
	bool do_it=
	  grid_x >= PBox[0] && grid_x < PBox[1] &&
	  grid_y >= PBox[2] && grid_y < PBox[3] &&
	  grid_z >= PBox[4] && grid_z < PBox[5];
	if(do_it && !period)
	  {
	    do_it=
	  grid_x >= moat && grid_x < glm &&
	  grid_y >= moat && grid_y < glm &&
	  grid_z >= moat && grid_z < glm;

	  }
	if(do_it)
	  {
	    int grid=fractal.where_1(grid_x,grid_y,grid_z);
	    assert(grid >= 0);
	    Point* p_point=group.list_points[grid];
	    p_point->list_particles.push_back(p);
	  }
	else
	  {
	    p->set_p_highest_level_group(0);
	    cout << " bad ins " << particle << endl;
	  }
      }
    cout << "end tree " << " " << group.list_points.size() << endl;
    if(fractal.get_level_max() > 0)
      group.set_number_high_groups(-1);
    else
      group.set_number_high_groups(0);
    cout << "exit treestart" << endl;
    t3=clock();
    cout << "tree time " << t1-t0 << " " << t2-t1 << " " << t3-t2 << endl;
  }
  //
}
