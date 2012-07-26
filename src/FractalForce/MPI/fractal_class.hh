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
    vector <bool>Periods;
    bool halo_fixed;
    unsigned int minimum_number;
    int level_max;
    int padding;
    //
    bool MPIrun;
    int TotalNodes;
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
    vector <int> Buffer;
    vector <int> Slice;
    vector <int> Box_to_Slices;
    vector <int> Box_from_Slices;
    vector <int> Slice_from_Boxes;
    vector <int> Slice_to_Boxes;
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
    int moat_0;
    int random_offset;
    int maxits;
    vector <clock_t> time_1;
    vector <clock_t> time_2;
    vector <clock_t> delta_time;
    vector <clock_t> total_time;
    vector <clock_t> time_g;
    vector <clock_t> delta_g;
    vector <clock_t> total_g;
    vector <clock_t> time_p;
    vector <clock_t> delta_p;
    vector <clock_t> total_p;
    bool debug;
    int highest_level_used;
    Fractal_Memory* p_generated_from;
    vector <string> time_string;
    int steps;
    double omega_start;
  public:
    Mess* p_mess;
    vector <Particle*> particle_list;
    vector <Particle*> particle_list_world;
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
    Fractal():
      number_particles(262144),
      number_particles_world(262144),
      grid_length(64),
      periodic(true),
      halo_fixed(false),
      minimum_number(8),
      level_max(8),
      padding(0),
      MPIrun(false),
      masks(0),
      epsilon_sor(6.0e-5),
      force_max(-1.0),
      halo_scale(1.0),
      halo_density0(1.0),
      density_0(0.0),
      moat_0(1),
      random_offset(0),
      maxits(250),
      debug(false),
      highest_level_used(0),
      omega_fraction(2.0/3.0)
    {
      p_mess=0;
      FractalNodes=get_FractalNodes();
      FractalRank=get_FractalRank();
      Box.assign(6,grid_length-1);
      Box[0]=0;
      Box[2]=0;
      Box[4]=0;
      //
      Periods.resize(3);
      Periods[0]=periodic && Box[0] == 0 && Box[1] == grid_length-1;
      Periods[1]=periodic && Box[2] == 0 && Box[3] == grid_length-1;
      Periods[2]=periodic && Box[4] == 0 && Box[5] == grid_length-1;
      Buffer.assign(6,0);
      BBox.assign(6,0);
      //      calc_Buffer();
      //      calc_BufferLev();
      //
      steps=0;
      omega_start=1.0;
      p_generated_from=0;
      time_1.assign(30,0);
      time_2.assign(30,0);
      delta_time.assign(30,0);
      total_time.assign(30,0);
      delta_g.assign(21,0);
      delta_p.assign(21,0);
      total_g.assign(21,0);
      total_p.assign(21,0);
      time_g.assign(21,0);
      time_p.assign(21,0);
      time_string.resize(30);
      time_string[0]="Initial isolated solver\t";
      time_string[1]="tree start\t";
      time_string[2]="edge trouble\t";
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
      time_string[21]="poisson solver\t";
      time_string[22]="force at point\t";
      time_string[23]="force at particle\t";
      time_string[24]=" ";
      time_string[25]="Halo Force Fixed";
      time_string[26]="clean up\t";
      time_string[27]="";
      time_string[28]="";
      time_string[29]="Everything\t";
      rad.assign(101,0.0);
      grow.assign(101,0.0);
    }
    template <class M> Fractal(M& mem):
      density_0(0.0),
      debug(false),
      highest_level_used(0),
      omega_fraction(2.0/3.0)
    {
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
      moat_0=mem.moat_0;
      random_offset=mem.random_offset;
      maxits=mem.maxits;
      debug=mem.debug;
      //
      p_mess=mem.p_mess;
      MPIrun=mem.MPIrun;
      TotalNodes=get_TotalNodes();
      FractalNodes=mem.FractalNodes;
      FractalNodes0=mem.FractalNodes0;
      FractalNodes1=mem.FractalNodes1;
      FractalNodes2=mem.FractalNodes2;
      FractalRank=get_FractalRank();
      assert(FractalRank<TotalNodes);
      Box=mem.Boxes[FractalRank];
      BBox=mem.BBoxes[FractalRank];
      PBox=mem.PBoxes[FractalRank];
      PBoxLength=mem.PBoxesLength[FractalRank];
      Buffer=mem.Buffers[FractalRank];
      BoxLev=mem.BoxesLev[FractalRank];
      BBoxLev=mem.BBoxesLev[FractalRank];
      Periods=mem.Periods[FractalRank];
      Slice=mem.p_mess->Slices[FractalRank];
      RealBox=mem.RealBoxes[FractalRank];
      //
      time_1.assign(30,0);
      time_2.assign(30,0);
      delta_time.assign(30,0);
      total_time.assign(30,0);
      delta_g.assign(21,0);
      delta_p.assign(21,0);
      total_g.assign(21,0);
      total_p.assign(21,0);
      time_g.assign(21,0);
      time_p.assign(21,0);
      time_string.resize(30);
      time_string[0]="Initial isolated solver\t";
      time_string[1]="tree start\t";
      time_string[2]="edge trouble\t";
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
      time_string[21]="poisson solver\t";
      time_string[22]="force at point\t";
      time_string[23]="force at particle\t";
      time_string[24]=" ";
      time_string[25]="Halo Force Fixed";
      time_string[26]="clean up\t";
      time_string[27]="";
      time_string[28]="";
      time_string[29]="Everything\t";
      masks=mem.masks;
      masks_level=mem.masks_level;
      masks_center_x=mem.masks_center_x;
      masks_center_y=mem.masks_center_y;
      masks_center_z=mem.masks_center_z;
      masks_rad_x=mem.masks_rad_x;
      masks_rad_y=mem.masks_rad_y;
      masks_rad_z=mem.masks_rad_z;
      masks_square=mem.masks_square;
      number_particles=mem.max_particles;;
      rad.assign(101,0.0);
      grow.assign(101,0.0);
      cout << "Making Fractal " << this << endl;
    }

    ~Fractal()
    {    
      cout << "Ending Fractal " << this << endl;
    };
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
    int get_TotalNodes()
    {
      return 1;
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
    void getPBoxLength(vector <int>& PBL)
    {
      PBL=PBoxLength;
    }
    void getBuffer(vector <int>& Bu)
    {
      Bu=Buffer;
    }
    void assign_buffer_edge_passive(Point& point,const int& level,bool& buff,bool& edge,bool& pass)
    {
      vector <int> pos(3);
      point.get_pos_point(pos);
      edge=false;
      buff=false;
      pass=
	pos[0]< BBoxLev[level][0] ||
	pos[0]> BBoxLev[level][1] ||
	pos[1]< BBoxLev[level][2] ||
	pos[1]> BBoxLev[level][3] ||
	pos[2]< BBoxLev[level][4] ||
	pos[2]> BBoxLev[level][5];
      if(pass)
	return;
      buff=
	pos[0]== BBoxLev[level][0] ||
	pos[0]== BBoxLev[level][1] ||
	pos[1]== BBoxLev[level][2] ||
	pos[1]== BBoxLev[level][3] ||
	pos[2]== BBoxLev[level][4] ||
	pos[2]== BBoxLev[level][5];
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
      pass=false;

      pass=
	n[0]==PBox[0] || n[0] == PBox[1] ||
	n[1]==PBox[2] || n[1] == PBox[3] ||
	n[2]==PBox[4] || n[2] == PBox[5];
      if(pass)
	return;
      buff=
	n[0]==BBox[0] || n[0] == BBox[1] ||
	n[1]==BBox[2] || n[1] == BBox[3] ||
	n[2]==BBox[4] || n[2] == BBox[5];
      if(buff)
	return;
      edge=
	n[0]==Box[0] || n[0] == Box[1] ||
	n[1]==Box[2] || n[1] == Box[3] ||
	n[2]==Box[4] || n[2] == Box[5];
      inside=
	(Periods[0] || (n[0] > BBox[0] && n[0] < BBox[1]))&&
	(Periods[1] || (n[1] > BBox[2] && n[1] < BBox[3]))&&
	(Periods[2] || (n[2] > BBox[4] && n[2] < BBox[5]));
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
      cout << "density_0 " << density_0 << endl;
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
    int get_moat_0()
    {
      return moat_0;
    }
    void set_moat_0(const int& m)
    {
      moat_0=m;
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
    void set_Periods(vector <bool>& Per)
    {   
      Periods=Per;
    }
    void get_Periods(vector <bool>& Per)
    {   
      Per=Periods;
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
    void timing_lev(const int& what,const int& level)
    {
      static ofstream FileTimeLev;
      if(!FileTimeLev.is_open())
	FileTimeLev.open("timing_lev.d");
      FileTimeLev.precision(2);
      if(what == -2)
	time_g[level]=clock();
      else if(what == -1)
	time_p[level]=clock();
      else if(what == 2)
	delta_g[level]=clock()-time_g[level];
      else if(what == 1)
	delta_p[level]=clock()-time_p[level];
      else if(what ==0)
	{
	  FileTimeLev.precision(2);
	  FileTimeLev << " " << endl;
	  FileTimeLev << " steps " << steps << endl;
	  for(int ni=0;ni<=level_max;ni++)
	    {
	      total_g[ni]+=delta_g[ni];
	      total_p[ni]+=delta_p[ni];
	      double dtg=delta_g[ni]/(double)CLOCKS_PER_SEC;
	      double dtp=delta_p[ni]/(double)CLOCKS_PER_SEC;
	      double totalg=total_g[ni]/(double)CLOCKS_PER_SEC;
	      double totalp=total_p[ni]/(double)CLOCKS_PER_SEC;
	      FileTimeLev << steps <<"\t" << ni << scientific << "\t" << dtg << "\t" << totalg << "\t" << dtp << "\t" << totalp << endl;
	    }
	}
      else
	assert(0);
    }
    void timing(const int& what, const int& which)
    {
      static ofstream FileTime;
      if(!FileTime.is_open())
	FileTime.open("timing.d");
      FileTime.precision(2);
      if(what == -1)
	time_1[which]=clock();
      else if(what == 1)
	{
	  time_2[which]=clock();
	  delta_time[which]+=time_2[which]-time_1[which];
	  if(which == 4)
	    delta_time[4]-=time_2[6]-time_1[6];
	  if(which == 1)
	    delta_time[1]-=time_2[2]-time_1[2];
	}
      else if(what == 0)
	{
	  steps++;
	  FileTime << " " << endl;
	  FileTime << " steps " << steps << endl;
	  for(int i=0; i < 30; i++)
	    total_time[i]+=delta_time[i];
	  double dt29=(delta_time[29]/(double)CLOCKS_PER_SEC)+1.0e-6;
	  double dtt29=(total_time[29]/(double)CLOCKS_PER_SEC)+1.0e-6;
	  for(int i=0; i < 30; i++)
	    {
	      double dt=delta_time[i]/(double)CLOCKS_PER_SEC;
	      double dtt=total_time[i]/(double)CLOCKS_PER_SEC;
	      FileTime << "timing " << steps << "\t" << i << " \t" << scientific << dt << "\t"  << dtt << "\t" ;
	      FileTime << fixed << 100.0*dt/dt29 << "\t" << 100.0*dtt/dtt29 << "\t" << time_string[i] << endl;
	    }
	}
      else if(what== -2)
	delta_time.assign(30,0);
    }
    void numbers3(Point* p_point,vector <int>& numbers)
    {
      int number=p_point->get_point_to_number();
      numbers[2]=number % PBoxLength[2];
      numbers[1]=(number/PBoxLength[2]) % PBoxLength[1];
      numbers[0]=number/(PBoxLength[1]*PBoxLength[2]);
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
      int kk=k;
      if(Periods[2])
	{
	  if(kk < PBox[4])
	    kk+=grid_length;
	  if(kk > PBox[5])
	    kk-=grid_length;
	}
      else
	if(kk < PBox[4] || kk > PBox[5]) return -1;
      int jj=j;
      if(Periods[1])
	{
	  if(jj < PBox[2])
	    jj+=grid_length;
	  if(jj > PBox[3])
	    jj-=grid_length;
	}
      else
	if(jj < PBox[2] || jj > PBox[3]) return -1;
      int ii=i;
      if(Periods[0])
	{
	  if(ii < PBox[0])
	    ii+=grid_length;
	  if(ii > PBox[1])
	    ii-=grid_length;
	}
      else
	if(ii < PBox[0] || ii > PBox[1]) return -1;
      return (kk-PBox[4])+((jj-PBox[2])+(ii-PBox[0])*PBoxLength[1])*PBoxLength[2];
    }
    int where(const int& nx,const int& ny,const int& nz,vector <int>& boX,vector <int>& boXL)
    {
      return (nz-boX[4])+((ny-boX[2])+(nx-boX[0])*boXL[1])*boXL[2];
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
	      pos[j]++;
	      doit=true;
	    }
	  else if(pos[j] >= 1.0)
	    {
	      pos[j]--;
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
