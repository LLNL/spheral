#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
//
namespace FractalSpace
{
  void tree_start_mini(Group& group, Fractal& fractal,Fractal_Memory& mem,Misc& misc)
  {
    fractal.timing(-1,1);
    ofstream& FileFractal=fractal.p_file->DUMPS;
    double a_grid_length=(double)fractal.get_grid_length();
    vector <int>Boxu(6);
    vector <Point*>ud(6);
    double t0,t1,t2,t3;
    t0=fractal.p_mess->Clock();
    FileFractal << "enter treestart" << "\n";
    int point_counter=0;
    mem.total_points_generated=0;
    mem.total_points_used=0;
    Point* p_point=0;
    Point* new_points=0;
    vector <double> RealBox;
    fractal.getRealBox(RealBox);
    Misc::vector_print(RealBox,FileFractal);
    vector <int> Box;
    fractal.getBox(Box);
    Misc::vector_print(Box,FileFractal);
    bool MPIrun=mem.MPIrun;
    group.set_buffer_group(MPIrun);
    vector <int> Buffer(6);
    fractal.getBuffer(Buffer);
    Misc::vector_print(Buffer,FileFractal);
    vector <int> BBox(6);
    fractal.getBBox(BBox);
    Misc::vector_print(BBox,FileFractal);
    vector <int> PBox(6);
    fractal.getPBox(PBox);
    Misc::vector_print(PBox,FileFractal);
    vector <int> PBoxLength(3);
    fractal.getPBoxLength(PBoxLength);
    Misc::vector_print(PBoxLength,FileFractal);
    int volume=PBoxLength[0]*PBoxLength[1]*PBoxLength[2];
    vector <bool>real_point(volume,false);
    vector <bool>it_is_a_point(volume,false);
    vector <Point*> pmyself(volume,Point::nothing);
    vector <double> pos(3);
    //
    for(int particle=0; particle < fractal.get_number_particles(); ++particle)
      {
	Particle* p=fractal.particle_list[particle];
	p->get_pos(pos);
	double a_nx=pos[0]*a_grid_length;
	double a_ny=pos[1]*a_grid_length;
	double a_nz=pos[2]*a_grid_length;
	int nx=a_nx;
	int ny=a_ny;
	int nz=a_nz;
	if(a_nx < 0.0)
	  nx=-static_cast<int>(-a_nx+1.0);
	if(a_ny < 0.0)
	  ny=-static_cast<int>(-a_ny+1.0);
	if(a_nz < 0.0)
	  nz=-static_cast<int>(-a_nz+1.0);
	int n=fractal.where_1(nx,ny,nz);
	assert(n>=0);
	real_point[n]=true;
      }
    it_is_a_point=real_point;
    for(int nz=PBox[4];nz <= PBox[5];nz++)
      {
	for(int ny=PBox[2];ny <= PBox[3];ny++)
	  {
	    for(int nx=PBox[0];nx <= PBox[1];nx++)
	      {
		if(real_point[fractal.where_1(nx,ny,nz)])
		  {
		    for(int dz=0;dz<=1;dz++)
		      {
			for(int dy=0;dy<=1;dy++)
			  {
			    for(int dx=0;dx<=1;dx++)
			      {
				int n=fractal.where_1(nx+dx,ny+dy,nz+dz);
				if(n>= 0) 
				  it_is_a_point[n]=true;
			      }
			  }
		      }
		  }
	      }
	  }
      }
    real_point=it_is_a_point;
    if(fractal.get_padding() != 0)
      {
	for(int nz=PBox[4];nz <= PBox[5];nz++)
	  {
	    for(int ny=PBox[2];ny <= PBox[3];ny++)
	      {
		for(int nx=PBox[0];nx <= PBox[1];nx++)
		  {
		    int n=fractal.where_1(nx,ny,nz);
		    if(real_point[n])
		      {
			for(int dz=-1;dz<=1;dz++)
			  {
			    for(int dy=-1;dy<=1;dy++)
			      {
				for(int dx=-1;dx<=1;dx++)
				  {
				    int n=fractal.where_1(nx+dx,ny+dy,nz+dz);
				    if(n >= 0)
				      it_is_a_point[n]=true;
				  }			      }
			  }
		      }
		  }
	      }
	  }
	real_point=it_is_a_point;
      }

    for(int nz=PBox[4];nz <= PBox[5];nz++)
      {
	for(int ny=PBox[2];ny <= PBox[3];ny++)
	  {
	    for(int nx=PBox[0];nx <= PBox[1];nx++)
	      {
		if(!real_point[fractal.where_1(nx,ny,nz)])
		  continue;
		int n=fractal.where_1(nx-1,ny,nz);
		if(n >= 0)
		  it_is_a_point[n]=true;
		n=fractal.where_1(nx+1,ny,nz);
		if(n >= 0)
		  it_is_a_point[n]=true;
		n=fractal.where_1(nx,ny-1,nz);
		if(n >= 0)
		  it_is_a_point[n]=true;
		n=fractal.where_1(nx,ny+1,nz);
		if(n >= 0)
		  it_is_a_point[n]=true;
		n=fractal.where_1(nx,ny,nz-1);
		if(n >= 0)
		  it_is_a_point[n]=true;
		n=fractal.where_1(nx,ny,nz+1);
		if(n >= 0)
		  it_is_a_point[n]=true;
	      }
	  }
      }
    int total_points=0;
    vector <int>total_pointsZ(PBox[5]+1,0);
    for(int nz=PBox[4];nz <= PBox[5];nz++)
      {
	for(int ny=PBox[2];ny <= PBox[3];ny++)
	  {
	    for(int nx=PBox[0];nx <= PBox[1];nx++)
	      {
		if(it_is_a_point[fractal.where_1(nx,ny,nz)])
		  {
		    total_pointsZ[nz]++;
		    total_points++;
		  }
	      }
	  }
      }
    //
    FileFractal << " To generate point OK in treestart mini " << mem.p_mess->FractalRank << " " << total_points << " " << volume << " " << fractal.get_number_particles() << "\n";
    group.set_group_number(0);
    vector <int>grid(3);
    bool inside,edge,buff,pass;
    int countit=0;
    for(int gridz=PBox[4];gridz <= PBox[5];gridz++)
      {
	int new_counter=0;
	try
	  {    
	    countit+=total_pointsZ[gridz];
	    new_points=new Point[total_pointsZ[gridz]];
	    group.list_new_points.push_back(new_points);
	    mem.total_points_generated=countit;
	  }
	catch(bad_alloc& ba)
	  {
	    cerr << " bad Point allocation in tree start mini " << mem.p_mess->FractalRank << " " << total_points << " " << gridz << " " << total_pointsZ[gridz] << " " << ba.what() << "\n";
	    cerr << " generated points bad in treestart mini " << mem.p_mess->FractalRank << " " << countit << " " << volume << " " << fractal.get_number_particles() << "\n";
	    cerr << " gen bad " << PBox[4] << " " << gridz << " " << PBox[5] << "\n";
	    int FR=mem.p_mess->FractalRank;
	    cerr << " crash res " << FR << " " << mem.PBoxes[FR][0] << " " << mem.PBoxes[FR][1] << " " << mem.PBoxes[FR][2] << " ";
	    cerr << mem.PBoxes[FR][3] << " " << mem.PBoxes[FR][4] << " " << mem.PBoxes[FR][5] << "\n";
	    for(int FR=0;FR < mem.p_mess->FractalNodes;FR++)
	      {
		cerr << " crash res " << FR << " " << mem.PBoxes[FR][0] << " " << mem.PBoxes[FR][1] << " " << mem.PBoxes[FR][2] << " ";
		cerr << mem.PBoxes[FR][3] << " " << mem.PBoxes[FR][4] << " " << mem.PBoxes[FR][5] << "\n";
	      }
	    cerr.flush();
	    assert(0);
	  }
	grid[2]=gridz;
	for(int gridy=PBox[2];gridy <= PBox[3];gridy++)
	  {
	    grid[1]=gridy;
	    for(int gridx=PBox[0];gridx <= PBox[1];gridx++)
	      {
		grid[0]=gridx;
		int n=fractal.where_1(gridx,gridy,gridz);
		if(!it_is_a_point[n])
		  continue;
		p_point=&new_points[new_counter];
		pmyself[n]=p_point;
		new_counter++;
		mem.total_points_used++;
		Point& point=*p_point;
		group.list_points.push_back(p_point);
		point.set_point_pointer(p_point);
		point.set_real_pointer(0);
		point.set_pos_point(gridx*misc.zoom,gridy*misc.zoom,gridz*misc.zoom);
		point.set_p_in_group(misc.p_group_0);
		point.set_p_in_high_group(0);
		fractal.inside_edge_buffer_pass(grid,inside,edge,buff,pass);
		inside=inside && real_point[n];
		point.set_inside(inside);
		point.set_edge_buffer_passive_point(edge,buff,pass);
		point.set_number_in_list(point_counter);
		point_counter++;
	      }
	  }
      }
    FileFractal << "aa " << point_counter << "\n";
    t1=fractal.p_mess->Clock();
    FileFractal << "treestart enter second part" << "\n";
    for(int grid_z=PBox[4];grid_z <= PBox[5];grid_z++)
      {
	for(int grid_y=PBox[2];grid_y <= PBox[3];grid_y++)
	  {
	    for(int grid_x=PBox[0];grid_x <= PBox[1];grid_x++)
	      {
		int n=fractal.where_1(grid_x,grid_y,grid_z);
		if(!it_is_a_point[n])
		  continue;
		ud.assign(6,Point::nothing);
		fractal.where_6(grid_x,grid_y,grid_z,Boxu);
		for(int ni=0;ni<6;ni++)
		  {
		    if(Boxu[ni] >= 0)
		      ud[ni]=pmyself[Boxu[ni]];
		  }
		pmyself[n]->set_point_ud(ud);
	      }
	  }
      }
    if(misc.get_debug())
      {
	check_for_edge_trouble(fractal);
	FileFractal << "checked for edge" << "\n";
      }
    t2=fractal.p_mess->Clock();
    FileFractal << "total number of particles " << fractal.get_number_particles() << "\n";
    int partsin=0;
    int partsout=0;
    for(int particle=0; particle < fractal.get_number_particles(); ++particle)
      {
	Particle* p=fractal.particle_list[particle];
	p->set_p_highest_level_group(misc.p_group_0);
	p->set_highest_level(0);
	p->get_pos(pos);
	//
	double a_grid_x=pos[0]*a_grid_length;
	double a_grid_y=pos[1]*a_grid_length;
	double a_grid_z=pos[2]*a_grid_length;
	int grid_x=a_grid_x;
	int grid_y=a_grid_y;
	int grid_z=a_grid_z;
	if(a_grid_x < 0.0)
	  grid_x=-static_cast<int>(-a_grid_x+1.0);
	if(a_grid_y < 0.0)
	  grid_y=-static_cast<int>(-a_grid_y+1.0);
	if(a_grid_z < 0.0)
	  grid_z=-static_cast<int>(-a_grid_z+1.0);
	bool do_it=
	  grid_x >= PBox[0] && grid_x < PBox[1] &&
	  grid_y >= PBox[2] && grid_y < PBox[3] &&
	  grid_z >= PBox[4] && grid_z < PBox[5];
	bool real_particle=do_it;
	if(do_it)
	  real_particle=
	    pos[0] >= RealBox[0] && pos[0] < RealBox[1] &&
	    pos[1] >= RealBox[2] && pos[1] < RealBox[3] &&
	    pos[2] >= RealBox[4] && pos[2] < RealBox[5];
	p->set_real_particle(real_particle);
	if(do_it)
	  {
	    int grid=fractal.where_1(grid_x,grid_y,grid_z);
	    assert(grid >= 0);
	    Point* p_point=pmyself[grid];
	    assert(p_point);
	    assert(real_point[grid]);
	    p_point->list_particles.push_back(p);
	    partsin++;
	  }
	else
	  {
	    p->set_p_highest_level_group(0);
	    partsout++;
	    FileFractal << " bad ins " << particle << "\t";
	    FileFractal << grid_x << "\t" << grid_y << "\t" << grid_z << "\t";
	    FileFractal << pos[0] << "\t" << pos[1] << "\t" << pos[2] << "\t"  << "\n";
	  }
      }
    FileFractal << "total number of particles in and out " << partsin << "\t" << partsout << "\n";
    FileFractal << "end tree " << "\t" << group.list_points.size() << "\n";
    if(fractal.get_level_max() > 0)
      group.set_number_high_groups(-1);
    else
      group.set_number_high_groups(0);
    FileFractal << "exit treestart" << "\n";
    t3=fractal.p_mess->Clock();
    FileFractal << "tree time " << t1-t0 << "\t" << t2-t1 << "\t" << t3-t2 << "\n";
    fractal.timing(1,1);
  }
}
