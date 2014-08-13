#ifndef _Fractal_Defined_
#define _Fractal_Defined_
namespace FractalSpace
{
  class Fractal 
  {
    int number_particles;
    int number_particles_world;
    int grid_length;
    bool periodic;
    bool halo_fixed;
    unsigned int minimum_number;
    int level_max;
    int padding;
    //
    bool MPIrun;
    int FractalNodes;
    int FractalNodes0;
    int FractalNodes1;
    int FractalNodes2;
    int FractalRank;
    vector <int> Box;
    vector <int> BBox;
    vector <int> PBox;
    vector <double> RealBox;
    vector <int> PBoxLength; 
    vector < vector <int> >BoxLev;
    vector < vector <int> >BBoxLev;
    vector < vector <int> >PBoxLev;
    vector <int> Buffer;
    double clocks_per_sec;
    //
    int masks;
    vector <double> masks_center_x;
    vector <double> masks_center_y;
    vector <double> masks_center_z;
    vector <double> masks_rad_x;
    vector <double> masks_rad_y;
    vector <double> masks_rad_z;
    vector <bool> masks_square;
    vector <int> masks_level;
    double epsilon_sor;
    double force_max;
    double halo_scale;
    double halo_density0;
    double density_0;
    int random_offset;
    int maxits;
    vector <double> time_1;
    vector <double> time_2;
    vector <double> delta_time;
    vector <double> total_time;
    vector <double> time_g;
    vector <double> delta_g;
    vector <double> total_g;
    vector <double> time_p;
    vector <double> delta_p;
    vector <double> total_p;
    bool debug;
    int highest_level_used;
    Fractal_Memory* p_generated_from;
    vector <string> time_string;
    int steps;
    double omega_start;
    double base_mass;
  public:
    Mess* p_mess;
    File* p_file;
    vector <Particle*> particle_list;
    vector <Particle*> particle_list_world;
    vector <Particle*> pseudo_particle_list;
    Particle* part_list_tmp;
    double omega_fraction;
    vector <double> rad;
    vector <double> grow;
    static bool first_time_solver;
    static string power_spec;
    static string integrator;
    static string energy_method;
    static string sim_parameters;
    static string force_fixed;
    static string particles;
    static string pot_solver;
    static string vel;
    Fractal()
    {
      //      cout << " made ghost fractal" << "\n";
    }
    template <class M> Fractal(M& mem):
      density_0(0.0),
      debug(false),
      highest_level_used(0),
      omega_fraction(2.0/3.0)
    {
      //      cout << " starting fractal " << "\n";
      //      clocks_per_sec=static_cast<double>(CLOCKS_PER_SEC);
      clocks_per_sec=1.0;
      steps=0;
      omega_start=mem.omega_start;
      p_generated_from=&mem;
      number_particles=mem.number_particles;
      number_particles_world=number_particles;
      grid_length=mem.grid_length;
      periodic=mem.periodic;
      halo_fixed=mem.halo_fixed;
      minimum_number=mem.minimum_number;
      level_max=mem.level_max;
      padding=mem.padding;
      epsilon_sor=mem.epsilon_sor;
      force_max=mem.force_max;
      halo_scale=mem.halo_scale;
      halo_density0=mem.halo_density0;
      random_offset=mem.random_offset;
      maxits=mem.maxits;
      debug=mem.debug;
      base_mass=mem.base_mass;
      //
      //      cout << " fractal start a " << FractalNodes0 << " " << FractalNodes1 << " " << FractalNodes2 << " " << FractalNodes << "\n";
      p_mess=mem.p_mess;
      p_file=mem.p_file;
      MPIrun=mem.MPIrun;
      FractalNodes=mem.FractalNodes;
      FractalNodes0=mem.FractalNodes0;
      FractalNodes1=mem.FractalNodes1;
      FractalNodes2=mem.FractalNodes2;
      //      cout << " fractal start b " << FractalNodes0 << " " << FractalNodes1 << " " << FractalNodes2 << " " << FractalNodes << "\n";
      //      cout << " fractal start c " << p_mess << " " << p_file << "\n";
      FractalRank=get_FractalRank();
      //      cout << FractalRank << "\n";
      assert(FractalRank<FractalNodes);
      Box=mem.Boxes[FractalRank];
      BBox=mem.BBoxes[FractalRank];
      PBox=mem.PBoxes[FractalRank];
      PBoxLength=mem.PBoxesLength[FractalRank];
      Buffer=mem.Buffers[FractalRank];
      BoxLev=mem.BoxesLev[FractalRank];
      BBoxLev=mem.BBoxesLev[FractalRank];
      PBoxLev=mem.PBoxesLev[FractalRank];
      RealBox=mem.RealBoxes[FractalRank];
      //
      p_file->FileFractal << "Box frac " << Box[0] << " " << Box[1] << " " << Box[2] << " " << Box[3] << " " << Box[4] << " " << Box[5] << "\n";
      //
      time_1.assign(50,0);
      time_2.assign(50,0);
      delta_time.assign(50,0);
      total_time.assign(50,0);
      delta_g.assign(21,0);
      delta_p.assign(21,0);
      total_g.assign(21,0);
      total_p.assign(21,0);
      time_g.assign(21,0);
      time_p.assign(21,0);
      time_string.resize(50);
      time_string[0]="Initial isolated solver\t";
      time_string[1]="tree start\t";
      time_string[2]="Edge Trouble\t";
      time_string[3]="assign density 0 ";
      time_string[4]="periodic solver\t";
      time_string[5]="isolated solver\t";
      time_string[6]="Power Spectrum\t";
      time_string[7]="Force at Point 0";
      time_string[8]="Force at Particle 0";
      time_string[9]="high points\t";
      time_string[10]="buffer points\t";
      time_string[11]="high pairs\t";
      time_string[12]="equivalence class\t";
      time_string[13]="high groups\t";
      time_string[14]="daughter group\t";
      time_string[15]="Connect Points\t";
      time_string[16]="Test Tree\t";
      time_string[17]="heavies\t";
      time_string[18]="particle lists\t";
      time_string[19]="assign density\t";
      time_string[20]="potential start\t";
      time_string[21]="\t";
      time_string[22]="force at point\t";
      time_string[23]="force at particle\t";
      time_string[24]="rho slice pot\t";
      time_string[25]="Halo Force Fixed";
      time_string[26]="clean up\t";
      time_string[27]="scatter particles\t";
      time_string[28]="gather particles\t";
      time_string[29]="\t";
      time_string[30]="\t";
      time_string[31]="Poisson Solver\t";
      time_string[32]="Hypre Search\t";
      time_string[33]="Wait Scatter\t";
      time_string[34]="Wait Dens-Slices\t";
      time_string[35]="Wait Slices-Pot\t";
      time_string[36]="Wait Hypre a\t";
      time_string[37]="Wait Hypre b\t";
      time_string[38]="Wait Gather\t";
      time_string[39]="Wait Info-Slices\t";
      time_string[40]="Wait Slices-PotIn\t";
      time_string[41]="Wait Global Level\t";
      time_string[42]="\t";
      time_string[43]="\t";
      time_string[44]="Tree Dump\t";
      time_string[45]="Cosmo Startup\t";
      time_string[46]="Tree Total\t";
      time_string[47]="Poisson Total\t";
      time_string[48]="Wait Time\t";
      time_string[49]="Everything\t";
      masks=mem.masks;
      masks_level=mem.masks_level;
      masks_center_x=mem.masks_center_x;
      masks_center_y=mem.masks_center_y;
      masks_center_z=mem.masks_center_z;
      masks_rad_x=mem.masks_rad_x;
      masks_rad_y=mem.masks_rad_y;
      masks_rad_z=mem.masks_rad_z;
      masks_square=mem.masks_square;
      rad.assign(101,0.0);
      grow.assign(101,0.0);
      //      cout << "Making Fractal " << this << "\n";
    }
    ~Fractal()
    {    
      //      cout << "Ending Fractal " << this << "\n";
    };
    void redo(Fractal_Memory* PFM)
    {
      Box=PFM->Boxes[FractalRank];
      BBox=PFM->BBoxes[FractalRank];
      PBox=PFM->PBoxes[FractalRank];
      PBoxLength=PFM->PBoxesLength[FractalRank];
      Buffer=PFM->Buffers[FractalRank];
      BoxLev=PFM->BoxesLev[FractalRank];
      BBoxLev=PFM->BBoxesLev[FractalRank];
      PBoxLev=PFM->PBoxesLev[FractalRank];
      RealBox=PFM->RealBoxes[FractalRank];
    }
    int get_FractalRank()
    {
      return p_mess->FractalRank;
    }
    int get_FractalNodes()
    {
      return p_mess->FractalNodes;
    }
    void set_omega_start(int& omega)
    {
      omega_start=omega;
    }
    double get_omega_start()
    {
      return omega_start;
    }
    double get_base_mass()
    {
      return base_mass;
    }
    void set_base_mass(double bm)
    {
      base_mass=bm;
    }
    void setBox(vector <int>& B)
    {
      Box=B;
    }
    void setBBox(vector <int>& BB)
    {
      BBox=BB;
    }
    void getPBox(vector <int>& PB)
    {
      PB=PBox;
    }
    void setPBox(vector <int>& PB)
    {
      PBox=PB;
      PBoxLength[0]=PBox[1]-PBox[0]+1;
      PBoxLength[1]=PBox[3]-PBox[2]+1;
      PBoxLength[2]=PBox[5]-PBox[4]+1;
  }
    void setBuffer(vector <int>& BB)
    {
      Buffer=BB;
    }
    void getBox(vector <int>& B)
    {
      B=Box;
    }
    void getBBox(vector <int>& BB)
    {
      BB=BBox;
    }
    void getBBoxLev(vector <int>& BB,const int& level)
    {
      BB=BBoxLev[level];
    }
    void getPBoxLev(vector <int>& PB,const int& level)
    {
      PB=PBoxLev[level];
    }
    void getPBoxLength(vector <int>& PBL)
    {
      PBL=PBoxLength;
    }
    void getBuffer(vector <int>& Bu)
    {
      Bu=Buffer;
    }
    void getRealBox(vector <double>& RB)
    {
      RB=RealBox;
    }
    void assign_edge_buffer_passive(Point& point,const int& level,bool& edge,bool& buff,bool& pass,bool& really)
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
    void setMPIrun(const bool& Mr)
    {
      MPIrun=Mr;
    }
    bool getMPIrun()
    {
      return MPIrun;
    }
    void inside_edge_buffer_pass(vector <int>& n,bool& inside,bool& edge,bool& buff,bool& pass)
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
    Fractal_Memory* get_p_generated_from()
    {
      return p_generated_from;
    }
    void set_p_generated_from (Fractal_Memory* p)
    {
      p_generated_from=p;
    }
    double get_halo_scale()
    {
      return halo_scale;
    }
    double get_halo_density0()
    {
      return halo_density0;
    }
    int get_grid_length()
    {
      return grid_length;
    }
    void set_grid_length(const int& i)
    {
      grid_length=i;
    }
    int get_highest_level_used()
    {
      return highest_level_used;
    }
    void set_highest_level_used(const int& i)
    {
      highest_level_used=i;
    }
    double get_density_0()
    {
      return density_0;
    }
    void set_density_0(const double& d)
    {
      density_0=d;
      p_file->FileFractal << "density_0 " << density_0 << "\n";
    }
    int get_number_particles()
    {
      return number_particles;
    }
    void set_number_particles(const int& i)
    {
      number_particles=i;
    }
    int get_number_particles_world()
    {
      return number_particles_world;
    }
    void set_number_particles_world(const int& i)
    {
      number_particles_world=i;
    }
    unsigned int get_minimum_number()
    {
      return minimum_number;
    }
    void set_minimum_number(const int& i)
    {
      minimum_number=i;
    }
    int get_level_max()
    {
      return level_max;
    }
    void set_level_max(const int& i)
    {
      level_max=i;
    }
    int get_padding()
    {
      return padding;
    }
    void set_padding(const int& i)
    {
      padding=i;
    }
    void set_epsilon_sor(const double& e)
    {
      epsilon_sor=e;
    }
    double get_epsilon_sor()
    {
      return epsilon_sor;
    }
    int get_random_offset()
    {
      return random_offset;
    }
    void set_maxits(const int& m)
    {
      maxits=m;
    }
    int get_maxits()
    {
      return maxits;
    }
    template <class M> void set_masks(M& mem)
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
    void set_pos_mask(const int& i,const double& x, const double& y, const double& z, const double& r1, const double& r2,const double& r3)
    {
      masks_center_x[i]=x;
      masks_center_y[i]=y;
      masks_center_z[i]=z;
      masks_rad_x[i]=r1;
      masks_rad_y[i]=r2;
      masks_rad_z[i]=r3;
    }
    double get_level_mask(const int& i)
    {
      return masks_level[i];
    }
    void set_level_mask(const int& i, const int& k)
    {
      masks_level[i]=k;
    }
    void get_mask(const int& m,vector <double>& cen,vector <double>& rad,bool& square)
    {
      cen[0]=masks_center_x[m];
      cen[1]=masks_center_y[m];
      cen[2]=masks_center_z[m];
      rad[0]=masks_rad_x[m];
      rad[1]=masks_rad_y[m];
      rad[2]=masks_rad_z[m];
      square=masks_square[m];
    }
    int get_number_masks()
    {
      return masks;
    }
    void set_number_masks(const int& i)
    {
      masks=i;
    }
    bool get_debug()
    {
      return debug;
    }
    void set_force_max(const double& f_max)
    {
      force_max=f_max;
    }
    double get_force_max()
    {
      return force_max;
    }
    void set_periodic(const bool& i)
    {
      periodic=i;
    }
    bool get_periodic()
    {
      return periodic;
    }
    bool get_halo_fixed()
    {
      return halo_fixed;
    }
    void set_steps(int& s)
    {
      steps=s;
    }
    void add_steps(int& s)
    {
      steps+=s;
    }
    int get_steps()
    {
      return steps;
    }
    void print_list(int what)
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
    void print_list_world(int what)
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
    double get_delta_time(int which)
    {
      return delta_time[which];
    }
    void get_total_times(vector <double>& TT)
    {
      TT=total_time;
    }
    void timing_lev(const int& what,const int& level)
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
    void timing(const int& what, const int& which)
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
    void where_6(const int& i,const int& j,const int& k,vector <int>& Boxu)
    {
      Boxu.assign(6,-1);
      Boxu[0]=where_1(i-1,j,k);
      Boxu[1]=where_1(i+1,j,k);
      Boxu[2]=where_1(i,j-1,k);
      Boxu[3]=where_1(i,j+1,k);
      Boxu[4]=where_1(i,j,k-1);
      Boxu[5]=where_1(i,j,k+1);
    }
    int where_1(const int& i,const int& j,const int& k)
    {
      if(k < PBox[4] || k > PBox[5]) return -1;
      if(j < PBox[2] || j > PBox[3]) return -1;
      if(i < PBox[0] || i > PBox[1]) return -1;
      return (i-PBox[0])+((j-PBox[2])+(k-PBox[4])*PBoxLength[1])*PBoxLength[0];
    }
    int where(const int& nx,const int& ny,const int& nz,vector <int>& boX,vector <int>& boXL)
    {
      return (nz-boX[4])+((ny-boX[2])+(nx-boX[0])*boXL[1])*(boXL[2]+2);
    }
    void wrap(const int& p)
    {
      wrap(particle_list[p]);
    }
    void wrap()
    {
      for(unsigned int p=0;p<particle_list.size();p++)
	wrap(p);
    }
    void wrap(Particle* p)
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
    static double my_rand(const double& rand_max)
    {
      return (double)(rand())/rand_max;
    }
    static double my_rand_not_zero(const double& rand_max)
    {
      return (double)(max(rand(),1))/rand_max;
    }
    template <class M> static bool equal(vector <M>& a,vector <M>& b) 
    {
      unsigned int sa=a.size();
      if(sa != b.size()); 
	return false;
      for(unsigned int i=0;i<sa;i++)
	{
	  if(a[i] == b[i])
	    continue;
	  return false;
	}
      return true;
    }
  };
}
#endif
