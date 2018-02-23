#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
//
namespace FractalSpace
{
  void Fractal::redo(Fractal_Memory* PFM)
  {
    Box=PFM->Boxes[FractalRank];
    BBox=PFM->BBoxes[FractalRank];
    PBox=PFM->PBoxes[FractalRank];
    PBoxLength=PFM->PBoxesLength[FractalRank];
    Buffer=PFM->Buffers[FractalRank];
    //      BoxLev=PFM->BoxesLev[FractalRank];
    //      BBoxLev=PFM->BBoxesLev[FractalRank];
    //      PBoxLev=PFM->PBoxesLev[FractalRank];
    BoxLev=PFM->FRBoxesLev;
    BBoxLev=PFM->FRBBoxesLev;
    PBoxLev=PFM->FRPBoxesLev;
    RealBox=PFM->RealBoxes[FractalRank];
  }
  int Fractal::get_FractalRank() const
  {
    return p_mess->FractalRank;
  }
  int Fractal::get_FractalNodes() const
  {
    return p_mess->FractalNodes;
  }
  void Fractal::set_omega_start(int& omega)
  {
    omega_start=omega;
  }
  double Fractal::get_omega_start() const
  {
    return omega_start;
  }
  double Fractal::get_base_mass() const
  {
    return base_mass;
  }
  void Fractal::set_base_mass(double bm)
  {
    base_mass=bm;
  }
  void Fractal::setBox(vector <int>& B)
  {
    Box=B;
  }
  void Fractal::setBBox(vector <int>& BB)
  {
    BBox=BB;
  }
  void Fractal::getPBox(vector <int>& PB) const
  {
    PB=PBox;
  }
  void Fractal::setPBox(vector <int>& PB)
  {
    PBox=PB;
    PBoxLength[0]=PBox[1]-PBox[0]+1;
    PBoxLength[1]=PBox[3]-PBox[2]+1;
    PBoxLength[2]=PBox[5]-PBox[4]+1;
  }
  void Fractal::setBuffer(vector <int>& BB)
  {
    Buffer=BB;
  }
  void Fractal::getBox(vector <int>& B) const
  {
    B=Box;
  }
  void Fractal::getBBox(vector <int>& BB) const
  {
    BB=BBox;
  }
  void Fractal::getBBoxLev(vector <int>& BB,const int& level) const
  {
    BB=BBoxLev[level];
  }
  void Fractal::getPBoxLev(vector <int>& PB,const int& level) const
  {
    PB=PBoxLev[level];
  }
  void Fractal::getPBoxLength(vector <int>& PBL) const
  {
    PBL=PBoxLength;
  }
  void Fractal::getBuffer(vector <int>& Bu) const
  {
    Bu=Buffer;
  }
  void Fractal::getRealBox(vector <double>& RB) const
  {
    RB=RealBox;
  }
  void Fractal::assign_edge_buffer_passive(Point& point,const int& level,bool& edge,bool& buff,bool& pass,bool& really) const
  {
    vector <int> pos(3);
    point.get_pos_point(pos);
    edge=false;
    buff=false;
    really=false;
    pass=
      pos[0]< BBoxLev[level][0] ||
      pos[0]> BBoxLev[level][1] ||
      pos[1]< BBoxLev[level][2] ||
      pos[1]> BBoxLev[level][3] ||
      pos[2]< BBoxLev[level][4] ||
      pos[2]> BBoxLev[level][5];
    if(pass)
      {
	really=
	  pos[0]< PBoxLev[level][0] ||
	  pos[0]> PBoxLev[level][1] ||
	  pos[1]< PBoxLev[level][2] ||
	  pos[1]> PBoxLev[level][3] ||
	  pos[2]< PBoxLev[level][4] ||
	  pos[2]> PBoxLev[level][5];
	return;
      }
    buff=
      pos[0] < BoxLev[level][0] ||
      pos[0] > BoxLev[level][1] ||
      pos[1] < BoxLev[level][2] ||
      pos[1] > BoxLev[level][3] ||
      pos[2] < BoxLev[level][4] ||
      pos[2] > BoxLev[level][5];
    if(buff)
      return;
    edge=
      pos[0]== BoxLev[level][0] ||
      pos[0]== BoxLev[level][1] ||
      pos[1]== BoxLev[level][2] ||
      pos[1]== BoxLev[level][3] ||
      pos[2]== BoxLev[level][4] ||
      pos[2]== BoxLev[level][5];
  }
  void Fractal::setMPIrun(const bool& Mr)
  {
    MPIrun=Mr;
  }
  bool Fractal::getMPIrun() const
  {
    return MPIrun;
  }
  void Fractal::inside_edge_buffer_pass(vector <int>& n,bool& inside,bool& edge,bool& buff,bool& pass) const
  {
    inside=false;
    edge=false;
    buff=false;
    pass=
      n[0] < BBox[0] || n[0] > BBox[1] ||
      n[1] < BBox[2] || n[1] > BBox[3] ||
      n[2] < BBox[4] || n[2] > BBox[5];
    //      if(pass)
    //	return;
    buff=
      n[0] < Box[0] || n[0] > Box[1] ||
      n[1] < Box[2] || n[1] > Box[3] ||
      n[2] < Box[4] || n[2] > Box[5];
    //      if(buff)
    //	return;
    edge=
      n[0]==Box[0] || n[0] == Box[1] ||
      n[1]==Box[2] || n[1] == Box[3] ||
      n[2]==Box[4] || n[2] == Box[5];
    inside=
      n[0] > BBox[0] && n[0] < BBox[1] &&
      n[1] > BBox[2] && n[1] < BBox[3] &&
      n[2] > BBox[4] && n[2] < BBox[5];
  }
  Fractal_Memory* Fractal::get_p_generated_from() const
  {
    return p_generated_from;
  }
  void Fractal::set_p_generated_from (Fractal_Memory* p)
  {
    p_generated_from=p;
  }
  double Fractal::get_halo_scale() const
  {
    return halo_scale;
  }
  double Fractal::get_halo_density0() const
  {
    return halo_density0;
  }
  int Fractal::get_grid_length() const
  {
    return grid_length;
  }
  void Fractal::set_grid_length(const int& i)
  {
    grid_length=i;
  }
  int Fractal::get_highest_level_used() const
  {
    return highest_level_used;
  }
  void Fractal::set_highest_level_used(const int& i)
  {
    highest_level_used=i;
  }
  double Fractal::get_density_0() const 
  {
    return density_0;
  }
  void Fractal::set_density_0(const double& d)
  {
    density_0=d;
    p_file->FileFractal << "density_0 " << density_0 << "\n";
  }
  int Fractal::get_number_particles() const
  {
    return number_particles;
  }
  void Fractal::set_number_particles(const int& i)
  {
    number_particles=i;
  }
  int Fractal::get_number_particles_world() const
  {
    return number_particles_world;
  }
  void Fractal::set_number_particles_world(const int& i)
  {
    number_particles_world=i;
  }
  unsigned int Fractal::get_minimum_number() const
  {
    return minimum_number;
  }
  void Fractal::set_minimum_number(const int& i)
  {
    minimum_number=i;
  }
  int Fractal::get_level_max() const
  {
    return level_max;
  }
  void Fractal::set_level_max(const int& i)
  {
    level_max=i;
  }
  int Fractal::get_padding() const
  {
    return padding;
  }
  void Fractal::set_padding(const int& i)
  {
    padding=i;
  }
  void Fractal::set_Hypre_TOL(const double& e)
  {
    HTOL=e;
  }
  double Fractal::get_Hypre_TOL() const
  {
    return HTOL;
  }
  int Fractal::get_random_offset() const
  {
    return random_offset;
  }
  void Fractal::set_maxits(const int& m)
  {
    maxits=m;
  }
  int Fractal::get_maxits() const
  {
    return maxits;
  }
  template <class M> void Fractal::set_masks(M& mem)
  {
    masks=mem.masks;
    masks_level=mem.masks_level;
    masks_center_x=mem.masks_center_x;
    masks_center_y=mem.masks_center_y;
    masks_center_z=mem.masks_center_z;
    masks_rad_x=mem.masks_rad_x;
    masks_rad_y=mem.masks_rad_y;
    masks_rad_z=mem.masks_rad_z;
    masks_square=mem.masks_square;
  }
  void Fractal::set_pos_mask(const int& i,const double& x, const double& y, const double& z, const double& r1, const double& r2,const double& r3)
  {
    masks_center_x[i]=x;
    masks_center_y[i]=y;
    masks_center_z[i]=z;
    masks_rad_x[i]=r1;
    masks_rad_y[i]=r2;
    masks_rad_z[i]=r3;
  }
  double Fractal::get_level_mask(const int& i) const
  {
    return masks_level[i];
  }
  void Fractal::set_level_mask(const int& i, const int& k)
  {
    masks_level[i]=k;
  }
  void Fractal::get_mask(const int& m,vector <double>& cen,vector <double>& rad,bool& square) const
  {
    cen[0]=masks_center_x[m];
    cen[1]=masks_center_y[m];
    cen[2]=masks_center_z[m];
    rad[0]=masks_rad_x[m];
    rad[1]=masks_rad_y[m];
    rad[2]=masks_rad_z[m];
    square=masks_square[m];
  }
  int Fractal::get_number_masks() const
  {
    return masks;
  }
  void Fractal::set_number_masks(const int& i)
  {
    masks=i;
  }
  bool Fractal::get_debug() const
  {
    return debug;
  }
  void Fractal::set_force_max(const double& f_max)
  {
    force_max=f_max;
  }
  double Fractal::get_force_max() const
  {
    return force_max;
  }
  void Fractal::set_periodic(const bool& i)
  {
    periodic=i;
  }
  bool Fractal::get_periodic() const
  {
    return periodic;
  }
  bool Fractal::get_halo_fixed() const
  {
    return halo_fixed;
  }
  void Fractal::set_steps(int& s)
  {
    steps=s;
  }
  void Fractal::add_steps(int& s)
  {
    steps+=s;
  }
  int Fractal::get_steps() const
  {
    return steps;
  }
  void Fractal::print_list(int what)
  {
    static int steppp=0;
    vector <double>pos(3);
    vector <double>pf(4);
    fprintf(p_file->PFPos," how many particles in list %d %d \n",particle_list.size(),what);
    for(unsigned int p=0;p<particle_list.size();p++)
      {
	particle_list[p]->get_pos(pos);
	if(what == 0)
	  fprintf(p_file->PFPos," PARTS%d %6d %10.6E %10.6E %10.6E \n",steppp,p,pos[0],pos[1],pos[2]);
	else
	  {
	    particle_list[p]->get_field_pf(pf);
	    fprintf(p_file->PFPos," PARTS%d %6d %10.6E %10.6E %10.6E %10.6E %10.6E %10.6E %10.6E \n",steppp,p,pos[0],pos[1],pos[2],
		    pf[0],pf[1],pf[2],pf[3]);
	  }
      }
    if(what == 0)
      steppp++;
  }
  void Fractal::print_list_world(int what)
  {
    static int steppp=0;
    vector <double>pos(3);
    vector <double>pf(4);
    fprintf(p_file->PFPos," how many particles in world list %d %d \n",particle_list_world.size(),what);
    for(unsigned int p=0;p<particle_list_world.size();p++)
      {
	particle_list_world[p]->get_pos(pos);
	if(what == 0)
	  fprintf(p_file->PFPos," WORLD%d %6d %10.6E %10.6E %10.6E \n",steppp,p,pos[0],pos[1],pos[2]);
	else
	  {
	    particle_list_world[p]->get_field_pf(pf);
	    fprintf(p_file->PFPos," WORLD%d %6d %10.6E %10.6E %10.6E %10.6E %10.6E %10.6E %10.6E \n",steppp,p,pos[0],pos[1],pos[2],
		    pf[0],pf[1],pf[2],pf[3]);
	  }
      }
    if(what == 0)
      steppp++;
  }
  double Fractal::get_delta_time(int which) const
  {
    return delta_time[which];
  }
  void Fractal::get_total_times(vector <double>& TT) const
  {
    TT=total_time;
  }
  void Fractal::timing_lev(const int& what,const int& level)
  {
    if(what == -2)
      time_g[level]=p_mess->Clock();
    else if(what == -1)
      time_p[level]=p_mess->Clock();
    else if(what == 2)
      delta_g[level]=p_mess->Clock()-time_g[level];
    else if(what == 1)
      delta_p[level]=p_mess->Clock()-time_p[level];
    else if(what ==0)
      {
	fprintf(p_file->PFTimeLev,"\n steps %5d \n",steps);
	for(int ni=0;ni<=level_max;ni++)
	  {
	    total_g[ni]+=delta_g[ni];
	    total_p[ni]+=delta_p[ni];
	    double dtg=delta_g[ni]/clocks_per_sec;
	    double dtp=delta_p[ni]/clocks_per_sec;
	    double totalg=total_g[ni]/clocks_per_sec;
	    double totalp=total_p[ni]/clocks_per_sec;
	    fprintf(p_file->PFTimeLev," %5d \t %3d \t %10.2E \t %10.2E \t %10.2E \t %10.2E \n",steps,ni,dtg,totalg,dtp,totalp);
	  }
	if(FractalRank == 0)
	  fflush(p_file->PFTimeLev);
      }
    else
      assert(0);
  }
  void Fractal::timing(const int& what, const int& which)
  {
    if(what == -1)
      time_1[which]=p_mess->Clock();
    else if(what == 1)
      {
	time_2[which]=p_mess->Clock();
	delta_time[which]+=time_2[which]-time_1[which];
      }
    else if(what == 0)
      {
	steps++;
	fprintf(p_file->PFTime,"\n steps %5d \n",steps);
	for(int i=0; i < 50; i++)
	  total_time[i]+=delta_time[i];
	double dt49=(delta_time[49]/clocks_per_sec)+1.0e-6;
	double dtt49=(total_time[49]/clocks_per_sec)+1.0e-6;
	total_time[48]=0.0;
	for(int ni=33;ni<=41;ni++)
	  {
	    delta_time[48]+=delta_time[ni];
	    total_time[48]+=total_time[ni];
	  }
	for(int i=0; i < 50; i++)
	  {
	    double dt=delta_time[i]/clocks_per_sec;
	    double dtt=total_time[i]/clocks_per_sec;
	    fprintf(p_file->PFTime,"timing %5d \t %3d \t %10.2E \t %10.2E \t",steps,i,dt,dtt);
	    fprintf(p_file->PFTime,"%10.2f \t %10.2f \t %s \n",100.0*dt/dt49,100.0*dtt/dtt49,time_string[i].c_str());
	  }
	if(FractalRank == 0)
	  fflush(p_file->PFTime);
      }
    else if(what== -2)
      delta_time.assign(50,0);
  }
  void Fractal::where_6(const int& i,const int& j,const int& k,vector <int>& Boxu) const
  {
    Boxu.assign(6,-1);
    Boxu[0]=where_1(i-1,j,k);
    Boxu[1]=where_1(i+1,j,k);
    Boxu[2]=where_1(i,j-1,k);
    Boxu[3]=where_1(i,j+1,k);
    Boxu[4]=where_1(i,j,k-1);
    Boxu[5]=where_1(i,j,k+1);
  }
  int Fractal::where_1(const int& i,const int& j,const int& k) const
  {
    if(k < PBox[4] || k > PBox[5]) return -1;
    if(j < PBox[2] || j > PBox[3]) return -1;
    if(i < PBox[0] || i > PBox[1]) return -1;
    return (i-PBox[0])+((j-PBox[2])+(k-PBox[4])*PBoxLength[1])*PBoxLength[0];
  }
  int Fractal::where(const int& nx,const int& ny,const int& nz,vector <int>& boX,vector <int>& boXL) const
  {
    return (nz-boX[4])+((ny-boX[2])+(nx-boX[0])*boXL[1])*(boXL[2]+2);
  }
  void Fractal::wrap(const int& p)
  {
    wrap(particle_list[p]);
  }
  void Fractal::wrap()
  {
    for(unsigned int p=0;p<particle_list.size();p++)
      wrap(p);
  }
  void Fractal::wrap(Particle* p)
  {
    vector <double> pos(3);
    p->get_pos(pos);
    bool doit=false;
    for(int j=0;j<3;j++)
      {
	if(pos[j] < 0.0)
	  {
	    pos[j]+=1.0;
	    doit=true;
	  }
	else if(pos[j] >= 1.0)
	  {
	    pos[j]-=1.0;
	    doit=true;
	  }
      }
    if(doit)
      p->set_pos(pos);
  }
}
