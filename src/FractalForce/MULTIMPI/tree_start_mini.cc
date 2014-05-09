#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
//
namespace FractalSpace
{
  void tree_start_mini(Group& group, Fractal& fractal,Fractal_Memory& mem,Misc& misc)
  {
    bool dumpit=mem.p_mess->FractalRank == mem.p_mess->FractalNodes-1;
    dumpit=false;
    int minimum=fractal.get_minimum_number();
    minimum=1;
    ofstream& FileFractal=fractal.p_file->DUMPS;
    ofstream* p_FilePoint=&fractal.p_file->DUMPS;
    double a_grid_length=(double)fractal.get_grid_length();
    vector <int>Boxu(6);
    vector <Point*>ud(6);
    double t0,t1,t2,t3;
    t0=fractal.p_mess->Clock();
    FileFractal << "enter treestart" << "\n";
    //    bool period=fractal.get_periodic();
    int point_counter=0;
    mem.total_points_generated=0;
    mem.total_points_used=0;
    Point* p_point=0;
    Point* new_points=0;
    //    FileFractal << " generating points in treestart" << "\n";
    //    FileFractal << " Boxes in node " << "\n";
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
    mini_Point* mP=new mini_Point[volume];
    if(dumpit)
      {
	cerr << " BOXESA " << Box[0] << " " << Box[1] << " " << Box[2] << " " << Box[3] << " " << Box[4] << " " << Box[5] << endl;
	cerr << " BOXESB " << BBox[0] << " " << BBox[1] << " " << BBox[2] << " " << BBox[3] << " " << BBox[4] << " " << BBox[5] << endl;
	cerr << " BOXESC " << PBox[0] << " " << PBox[1] << " " << PBox[2] << " " << PBox[3] << " " << PBox[4] << " " << PBox[5] << endl;
	cerr << " VOLUMEA " << volume << endl;
      }
    for(int nz=PBox[4];nz <= PBox[5];nz++)
      for(int ny=PBox[2];ny <= PBox[3];ny++)
	for(int nx=PBox[0];nx <= PBox[1];nx++)
	  {
	    int n=fractal.where_1(nx,ny,nz);
	    mini_Point* pmP=&mP[n];
	    pmP->realpoint=false;
	    pmP->it_is_a_point=false;
	    pmP->numbers=0;
	    pmP->pmyself=0;
	  }
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
	(&mP[n])->realpoint=true;
	(&mP[n])->it_is_a_point=true;
	(&mP[n])->numbers++;
      }
    int na=0;
    for(int nz=PBox[4];nz <= PBox[5];nz++)
      {
	for(int ny=PBox[2];ny <= PBox[3];ny++)
	  {
	    for(int nx=PBox[0];nx <= PBox[1];nx++)
	      {
		if(!(&mP[fractal.where_1(nx,ny,nz)])->it_is_a_point)
		  continue;
		na++;
		if(dumpit)
		  cerr << " dumpA " << nx << " " << ny << " " << nz << " " << fractal.where_1(nx,ny,nz) << " " << na << endl;
	      }
	  }
      }
    //
    for(int nz=PBox[4];nz <= PBox[5];nz++)
      {
	for(int ny=PBox[2];ny <= PBox[3];ny++)
	  {
	    for(int nx=PBox[0];nx <= PBox[1];nx++)
	      {
		if((&mP[fractal.where_1(nx,ny,nz)])->realpoint)
		  {
		    for(int dz=0;dz<=1;dz++)
		      {
			for(int dy=0;dy<=1;dy++)
			  {
			    for(int dx=0;dx<=1;dx++)
			      {
				int n=fractal.where_1(nx+dx,ny+dy,nz+dz);
				if(n < 0) 
				  continue;
				(&mP[n])->it_is_a_point=true;
			      }
			  }
		      }
		  }
	      }
	  }
      }
    na=0;
    for(int nz=PBox[4];nz <= PBox[5];nz++)
      {
	for(int ny=PBox[2];ny <= PBox[3];ny++)
	  {
	    for(int nx=PBox[0];nx <= PBox[1];nx++)
	      {
		int n=fractal.where_1(nx,ny,nz);
		if(!(&mP[n])->it_is_a_point)
		  continue;
		(&mP[n])->realpoint=true;
		na++;
		if(dumpit)
		  cerr << " dumpB " << nx << " " << ny << " " << nz << " " << n << " " << na << endl;
	      }
	  }
      }
    if(fractal.get_padding() != 0)
      {
	for(int nz=PBox[4];nz <= PBox[5];nz++)
	  {
	    for(int ny=PBox[2];ny <= PBox[3];ny++)
	      {
		for(int nx=PBox[0];nx <= PBox[1];nx++)
		  {
		    int n=fractal.where_1(nx,ny,nz);
		    if((&mP[n])->realpoint)
		      {
			for(int dz=-1;dz<=1;dz++)
			  {
			    for(int dy=-1;dy<=1;dy++)
			      {
				for(int dx=-1;dx<=1;dx++)
				  {
				    int n=fractal.where_1(nx+dx,ny+dy,nz+dz);
				    if(n < 0)
				      continue;
				    (&mP[n])->it_is_a_point=true;
				  }			      }
			  }
		      }
		  }
	      }
	  }
	na=0;
	for(int nz=PBox[4];nz <= PBox[5];nz++)
	  {
	    for(int ny=PBox[2];ny <= PBox[3];ny++)
	      {
		for(int nx=PBox[0];nx <= PBox[1];nx++)
		  {
		    int n=fractal.where_1(nx,ny,nz);
		    if(!(&mP[n])->it_is_a_point)
		      continue;
		    (&mP[n])->realpoint=true;
		    na++;
		    if(dumpit)
		      cerr << " dumpC " << nx << " " << ny << " " << nz << " " << fractal.where_1(nx,ny,nz) << " " << na << endl;
		  }
	      }
	  }
      }

    for(int nz=PBox[4];nz <= PBox[5];nz++)
      {
	for(int ny=PBox[2];ny <= PBox[3];ny++)
	  {
	    for(int nx=PBox[0];nx <= PBox[1];nx++)
	      {
		if(!(&mP[fractal.where_1(nx,ny,nz)])->realpoint)
		  continue;
		int n=fractal.where_1(nx-1,ny,nz);
		if(n >= 0)
		  (&mP[n])->it_is_a_point=true;
		n=fractal.where_1(nx+1,ny,nz);
		if(n >= 0)
		  (&mP[n])->it_is_a_point=true;
		n=fractal.where_1(nx,ny-1,nz);
		if(n >= 0)
		  (&mP[n])->it_is_a_point=true;
		n=fractal.where_1(nx,ny+1,nz);
		if(n >= 0)
		  (&mP[n])->it_is_a_point=true;
		n=fractal.where_1(nx,ny,nz-1);
		if(n >= 0)
		  (&mP[n])->it_is_a_point=true;
		n=fractal.where_1(nx,ny,nz+1);
		if(n >= 0)
		  (&mP[n])->it_is_a_point=true;
	      }
	  }
      }
    na=0;
    for(int nz=PBox[4];nz <= PBox[5];nz++)
      {
	for(int ny=PBox[2];ny <= PBox[3];ny++)
	  {
	    for(int nx=PBox[0];nx <= PBox[1];nx++)
	      {
		if(!(&mP[fractal.where_1(nx,ny,nz)])->it_is_a_point)
		  continue;
		na++;
		if(dumpit)
		  cerr << " dumpD " << nx << " " << ny << " " << nz << " " << fractal.where_1(nx,ny,nz) << " " << na << endl;
	      }
	  }
      }
    int total_points=0;
    for(int nz=PBox[4];nz <= PBox[5];nz++)
      {
	for(int ny=PBox[2];ny <= PBox[3];ny++)
	  {
	    for(int nx=PBox[0];nx <= PBox[1];nx++)
	      {
		if((&mP[fractal.where_1(nx,ny,nz)])->it_is_a_point)
		  {
		    total_points++;
		  }
	      }
	  }
      }
    //
    cout << " To generate point OK in treestart mini " << mem.p_mess->FractalRank << " " << total_points << " " << volume << " " << fractal.get_number_particles() << "\n";
    try
      {    
	new_points=new Point[total_points];
      }
    catch(bad_alloc& ba)
      {
	cerr << " bad allocation in tree start mini " << mem.p_mess->FractalRank;
	cerr << " " << total_points << " ";
	cerr << ba.what() << "\n";
	mem.p_file->FlushAll();
	cerr << " generated points bad in treestart mini " << mem.p_mess->FractalRank << " " << total_points << " " << volume << " " << fractal.get_number_particles() << "\n";
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
    mem.total_points_generated=total_points;
    group.set_group_number(0);
    group.list_new_points.push_back(new_points);
    FileFractal << " generated points in treestart " << mem.p_mess->FractalRank << " " << total_points << " " << volume << " " << fractal.get_number_particles() << "\n";
    //
    //
    vector <int>grid(3);
    bool inside,edge,buff,pass;
    int new_counter=0;
    for(int gridz=PBox[4];gridz <= PBox[5];gridz++)
      {
	grid[2]=gridz;
	for(int gridy=PBox[2];gridy <= PBox[3];gridy++)
	  {
	    grid[1]=gridy;
	    for(int gridx=PBox[0];gridx <= PBox[1];gridx++)
	      {
		grid[0]=gridx;
		int n=fractal.where_1(gridx,gridy,gridz);
		if(!(&mP[n])->it_is_a_point)
		  continue;
		p_point=&new_points[new_counter];
		(&mP[n])->pmyself=p_point;
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
		inside=inside && (&mP[n])->realpoint;
		point.set_inside(inside);
		point.set_edge_buffer_passive_point(edge,buff,pass);
		point.set_number_in_list(point_counter);
		point.set_FILE(p_FilePoint);
		//		FileFractal << " insides " << gridx << " " << gridy << " " << gridz << " " << inside << "\n";
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
		if(!(&mP[n])->it_is_a_point)
		  continue;
		ud.assign(6,Point::nothing);
		fractal.where_6(grid_x,grid_y,grid_z,Boxu);
		for(int ni=0;ni<6;ni++)
		  {
		    if(Boxu[ni] < 0)
		      continue;
		    ud[ni]=(&mP[Boxu[ni]])->pmyself;
		  }
		(&mP[n])->pmyself->set_point_ud(ud);
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
	//
	//	int grid_x=static_cast<int>(floor(pos[0]*a_grid_length));
	//	int grid_y=static_cast<int>(floor(pos[1]*a_grid_length));
	//	int grid_z=static_cast<int>(floor(pos[2]*a_grid_length));
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
	    Point* p_point=(&mP[grid])->pmyself;
	    assert(p_point);
	    assert((&mP[grid])->realpoint);
	    p_point->list_particles.push_back(p);
	    partsin++;
	    //	    FileFractal << " doit " << particle << " " << p << " " << p_point << " " << grid << " " << grid_x << " " << grid_y << " " << grid_z << "\n";
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
    delete [] mP;
    mP=0;
    FileFractal << "total number of particles in and out " << partsin << "\t" << partsout << "\n";
    FileFractal << "end tree " << "\t" << group.list_points.size() << "\n";
    if(fractal.get_level_max() > 0)
      group.set_number_high_groups(-1);
    else
      group.set_number_high_groups(0);
    FileFractal << "exit treestart" << "\n";
    t3=fractal.p_mess->Clock();
    FileFractal << "tree time " << t1-t0 << "\t" << t2-t1 << "\t" << t3-t2 << "\n";
  }
  //
}
