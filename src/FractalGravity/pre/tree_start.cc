#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
//
namespace FractalSpace
{
  void tree_start(Group& group, Fractal& fractal,Fractal_Memory& mem,Misc& misc)
  {
    bool halo=fractal.get_halo();
    vector <Point*>ud(6);
    double t0,t1,t2,t3,t4;
    t0=clock();
    cout << "enter treestart" << endl;
    bool period=fractal.get_periodic();
    int length=fractal.get_grid_length();
    int length_ratio=fractal.get_length_ratio();
    int point_counter=0;
    //
    if(!period)
      length++;
    //
    int lengthm1=length-1;
    int length_2=length*length;
    int length_3=length_2*length;
    int length_z=length;
    if(length_ratio != 1 && !period) 
      {
	length_z=(length-1)/length_ratio+1;
	length_3=length_2*length_z;
      }
    int length_zm1=length_z-1;
    double a_grid_length=(double)fractal.get_grid_length();
    if(!mem.remember_points)
      mem.total_points_generated=0;
    mem.total_points_used=0;
    Point* p_point=0;
    static Point* new_points=0;
    int new_counter=0;
    bool ins_x=period;
    //
    if(mem.total_points_generated == 0)
      {
	cout << " generating points in treestart" << endl;
	new_points=new (nothrow) Point[length_3];
	assert(new_points);
	mem.total_points_generated=length_3;
	group.list_new_points.push_back(new_points);
	group.list_points.reserve(length_3);
	cout << " generated points in treestart" << endl;
      }
    //
    for(int grid=0;grid<length_3;grid++)
      {
	//	cout << "tree grid " << grid << endl;
	int grid_z=grid/length_2;
	int grid_y=(grid/length) % length;
	int grid_x=grid % length;
	if(!period)
	  {
	    ins_x=grid_z > 0 && grid_z < length_zm1 && grid_y > 0 && grid_y < lengthm1
	      && grid_x > 0 && grid_x < lengthm1;
	  }
	p_point=&new_points[new_counter];
	//	cout << "tree grid a " << grid << " " << p_point << " " << new_counter << endl;
	new_counter++;
	mem.total_points_used++;
	Point& point=*p_point;
	group.list_points.push_back(p_point);
	//	cout << "tree grid b " << grid << " "  << point_counter << " " << p_point << endl;
	point.set_point_to_number(point_counter);
	point_counter++;
	point.set_point_pointer(p_point);
	point.set_real_pointer(0);
	//	cout << "tree grid c " << grid << endl;
	point.set_pos_point(grid_x*misc.zoom,grid_y*misc.zoom,grid_z*misc.zoom);
	point.set_inside(ins_x);
	if(!ins_x) point.set_real_pointer(-1);
	//	cout << "tree grid d " << grid << endl;
	point.set_p_in_group(misc.p_group_0);
	point.set_p_in_high_group(0);
	//	cout << "tree grid e " << grid << endl;
      }
    t1=clock();
    cout << "treestart enter second part" << endl;
    int count=0;
    for(int grid=0;grid<length_3;grid++)
      {
	int grid_z=grid/length_2;
	int grid_y=(grid/length) % length;
	int grid_x=grid % length;
	Point& point=*group.list_points[grid];
	int up_x=where_3(grid_x+1,grid_y,grid_z,length,period);
	int up_y=where_3(grid_x,grid_y+1,grid_z,length,period);
	int up_z=where_3(grid_x,grid_y,grid_z+1,length,period);
	int down_x=where_3(grid_x-1,grid_y,grid_z,length,period);
	int down_y=where_3(grid_x,grid_y-1,grid_z,length,period);
	int down_z=where_3(grid_x,grid_y,grid_z-1,length,period);
	if(period || point.get_inside())
	  {
	    ud[0]=group.list_points[down_x];
	    ud[1]=group.list_points[up_x];
	    ud[2]=group.list_points[down_y];
	    ud[3]=group.list_points[up_y];
	    ud[4]=group.list_points[down_z];
	    ud[5]=group.list_points[up_z];
	    point.set_point_ud(ud);
	  }
	else
	  {
	    if(grid_x == lengthm1)
	      point.set_point_up_x(0);
	    else
	      point.set_point_up_x(ud[1]);
	    if(grid_y == lengthm1)
	      point.set_point_up_y(0);
	    else
	      point.set_point_up_y(ud[3]);
	    if(grid_z == length_zm1)
	      point.set_point_up_z(0);
	    else
	      point.set_point_up_z(ud[5]);
	    if(grid_x  == 0)
	      point.set_point_down_x(0);
	    else
	      point.set_point_down_x(ud[0]);
	    if(grid_y == 0)
	      point.set_point_down_y(0);
	    else
	      point.set_point_down_y(ud[2]);
	    if(grid_z == 0)
	      point.set_point_down_z(0);
	    else
	      point.set_point_down_z(ud[4]);
	  }
	count++;
      }
    t2=clock();
    cout << "treestart exit second part" << endl;
    if(misc.get_debug())
      check_for_edge_trouble(fractal);
    t3=clock();
    //
    cout << "checked for edge" << endl;
    int g_l=fractal.get_grid_length();
    int g_l_z=g_l/length_ratio;
    int moat=0;
    if(!period) moat=max(max(1,fractal.get_moat_0()),fractal.get_padding());
    double radmin2=0.0;
    if(halo) radmin2=pow(0.5-(int)moat/(int)fractal.get_grid_length(),2);
    vector <double> pos(3);
    bool do_it;
    int grid_x,grid_y,grid_z;
    for(int particle=0; particle < fractal.get_number_particles(); ++particle)
      {
	Particle* p=fractal.particle_list[particle];
	p->set_p_highest_level_group(misc.p_group_0);
	p->set_highest_level(0);
	do_it=period;
	if(period) wrap(*p);
	p->get_pos(pos);
	grid_x=(int)floor(pos[0]*a_grid_length);
	grid_y=(int)floor(pos[1]*a_grid_length);
	grid_z=(int)floor(pos[2]*a_grid_length);
	if(period)
	  {
	    assert(grid_x >= 0 && grid_x < length);
	    assert(grid_y >= 0 && grid_y < length);
	    assert(grid_z >= 0 && grid_z < length);
	  }
	else
	  {
	    if(halo)
	      do_it=p->get_r2(0.5,0.5,0.5) < radmin2;
	    else
	      {
		do_it=grid_x >= moat && grid_y >= moat && grid_z >=  moat &&
		  grid_x < g_l-moat && grid_y < g_l-moat && grid_z < g_l_z-moat;
	      }
	  }
	if(do_it)
	  {
	    int grid=where_3(grid_x,grid_y,grid_z,length,period);
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
    t4=clock();
    cout << "tree time " << t1-t0 << " " << t2-t1 << " " << t3-t2 << " " << t4-t3 << " " << endl;
  }
  //
  inline int where_3(const int& i, const int& j, const int& k, const int& m, const bool& periodic)
  {
    if(periodic)
      return (i+m) % m + ((j+m) % m + ((k+m) % m)*m)*m;
    return i+(j+k*m)*m;
  }
  void wrap(double& x,double& y,double& z)
  {
    wrap(x);
    wrap(y);
    wrap(z);
  }
  void wrap(int& x,int& y,int& z,int& width)
  {
    wrap(x,width);
    wrap(y,width);
    wrap(z,width);
  }
  void wrap(int& x,int& width)
  {
    if(x < 0)
      x+=width;
    else if(x >= width)
      x=x% width;
  }
  void wrap(double& x)
  {
    if(x < 0.0)
      x+=1.0;
    else if(x >= 1.0)
      x-=1.0;
  }
  void wrap(Particle& part)
  {
    double x,y,z;
    part.get_pos(x,y,z);
    wrap(x,y,z);
    part.set_pos(x,y,z);
  }
}
