#ifndef _Fractal_Memory_Defined_
#define _Fractal_Memory_Defined_
namespace FractalSpace
{
  class Fractal_Memory
  {
  public:
    //
    string BaseDirectory;
    string RUN;
    MPI_Comm FractalWorld;
    bool standalone;
    bool MPIrun;
    bool balance;
    int FractalNodes;
    int FractalNodes0;
    int FractalNodes1;
    int FractalNodes2;
    int FFTNodes;
    bool time_trial;
    int min_hypre_group_size;
    vector < vector <int> > Boxes;
    vector < vector <int> > BBoxes;
    vector < vector <int> > PBoxes;
    vector < vector <int> > Buffers;
    vector < vector < vector <int> > > BoxesLev;
    vector < vector < vector <int> > > BBoxesLev;
    vector < vector < vector <int> > > PBoxesLev;
    vector < vector < vector <int> > > HRBoxesLev;
    vector < vector < vector <int> > > HSBoxesLev;
    vector < vector <int> > PBoxesLength;
    vector < vector <double> > RealBoxes;
    vector < vector <double> > RealPBoxes;
    vector <int>ij_offsets;
    vector <int>ij_counts;
    string hypre_solver;
    string hypre_precond;
    int global_level_max;
    //    vector <double>total_time;
    //
    bool split_particles;
    bool amnesia;
    bool mind_wipe;
    bool fixed_potential;
    bool calc_shear;
    bool start_up;
    bool calc_density_particle;
    bool do_vel;
    bool do_var; 
    bool periodic;
    bool random_initial;
    bool debug;
    bool halo_fixed;
    bool momentum_conserve;
    int total_points_counter;
    int total_points_used;
    int total_points_generated;
    int new_points_gen;
    int random_gen;
    int highest_level_init;
    int norm_what;
    int spectrum_number;
    int number_particles;
    int grid_length;
    int moat_0;
    unsigned int minimum_number;
    int padding;
    int level_max;
    int number_steps_total;
    int number_steps_out;
    int random_offset;
    int maxits;
    double base_mass;
    double epsilon_sor;
    double force_max;
    double halo_scale;
    double halo_density0;
    double off_x;
    double off_y;
    double off_z;
    double sigma_initial;
    double scaling;
    double norm_scale;
    double power_slope;
    double cut_off;
    double cut_off_init;
    double pexp;
    double step_length;
    double omega_start;
    double lambda_start;
    double redshift_start;
    double omega_0;
    double omega_lambda;
    double omega_b;
    double box_length;
    double h;
    double sigma_0;
    //
    double udda;
    double potential_energy;
    double potential_energy_old;
    double kinetic_energy;
    double kinetic_energy_old;
    double arad;
    double time;
    double total_mass;
    int steps;
    //
    int crash_levels;
    double crash_pow;
    double density_crash;
    int max_particles;
    int splits;
    vector <double> splits_center_x;
    vector <double> splits_center_y;
    vector <double> splits_center_z;
    vector <double> splits_rad_x;
    vector <double> splits_rad_y;
    vector <double> splits_rad_z;
    vector <bool> splits_square;
    vector <int> splits_many;
    int masks;
    vector <double> masks_center_x;
    vector <double> masks_center_y;
    vector <double> masks_center_z;
    vector <double> masks_rad_x;
    vector <double> masks_rad_y;
    vector <double> masks_rad_z;
    vector <bool> masks_square;
    vector <int> masks_level;
    int masks_init;
    vector <double> masks_center_x_init;
    vector <double> masks_center_y_init;
    vector <double> masks_center_z_init;
    vector <double> masks_rad_x_init;
    vector <double> masks_rad_y_init;
    vector <double> masks_rad_z_init;
    vector <bool> masks_square_init;
    vector <int> masks_level_init;
    //
    vector < vector<Group*> > all_groups;
    Misc* p_misc; 
    Fractal* p_fractal;
    Mess* p_mess;
    File* p_file;
    //
    Fractal_Memory():
      //
      // default values
      // replace with your own values
      //
      BaseDirectory("FFRRAACCTTAALL/"),
      RUN("abc"),
      FractalWorld(MPI_COMM_WORLD),
      standalone(true),
      MPIrun(false),
      balance(0),
      FractalNodes(1),
      FractalNodes0(1),
      FractalNodes1(1),
      FractalNodes2(1),
      FFTNodes(1234567),
      time_trial(true),
      min_hypre_group_size(1),
      amnesia(true),
      mind_wipe(false),
      fixed_potential(false),
      calc_shear(false),
      start_up(false),
      calc_density_particle(false),
      do_vel(false),
      do_var(false), 
      periodic(false),
      random_initial(false),
      debug(false),
      halo_fixed(false),
      momentum_conserve(true),
      total_points_counter(0),
      total_points_used(0),
      total_points_generated(0),
      new_points_gen(9),
      random_gen(12345),
      highest_level_init(6),
      norm_what(0),
      spectrum_number(0),
      number_particles(262144),
      grid_length(128),
      moat_0(1),
      minimum_number(8),
      padding(-1),
      level_max(8),
      number_steps_total(113),
      number_steps_out(20),
      random_offset(0),
      maxits(20),
      base_mass(1.0),
      epsilon_sor(1.0e-7),
      force_max(-1.0),
      halo_scale(1.0),
      halo_density0(1.0),
      off_x(0.5),
      off_y(0.5),
      off_z(0.5),
      sigma_initial(-1.0),
      scaling(1.0),
      norm_scale(0.2),
      power_slope(-1.0),
      cut_off(1.0e6),
      cut_off_init(16.0),
      pexp(0.77),
      step_length(0.01),
      omega_start(-1.0),
      lambda_start(-1.0),
      redshift_start(49.0),
      omega_0(0.3),
      omega_lambda(0.7),
      omega_b(0.03),
      box_length(100.0),
      h(0.7),
      sigma_0(1.0),
    //
      udda(0.0),
      potential_energy(0.0),
      potential_energy_old(0.0),
      kinetic_energy(0.0),
      kinetic_energy_old(0.0),
      arad(1.0),
      time(2.0/3.0),
      total_mass(1.0),
      steps(0),
    //
      crash_levels(0),
      crash_pow(1.0),
      density_crash(1.0e30),
      max_particles(0),
      splits(0),
      masks(0)
    {
      p_misc=0;
      p_fractal=0;
      p_file=0;
      hypre_solver="AMG";
      hypre_precond="AMG";
      global_level_max=level_max;
      padding=min(padding,1);
      split_particles= force_max > 0.0;
      //      total_time.assign(50,0.0);
      //
    }
    ~Fractal_Memory()
    {
      cout << "Ending Fractal_Memory " << this << endl;
    }
    void calc_FractalNodes()
    {
      FractalNodes=FractalNodes0*FractalNodes1*FractalNodes2;
      MPIrun=FractalNodes > 1;
      Boxes.resize(FractalNodes);
      int length=grid_length;
      if(!periodic)
	length++;
      int count=0;
      int j2b=-1;
      for(int m2=0;m2<FractalNodes2;m2++)
	{
	  int j2a=j2b+1;
	  j2b=((m2+1)*length)/FractalNodes2-1;
	  int j1b=-1;
	  for(int m1=0;m1<FractalNodes1;m1++)
	    {
	      int j1a=j1b+1;
	      j1b=((m1+1)*length)/FractalNodes1-1;
	      int j0b=-1;
	      for(int m0=0;m0<FractalNodes0;m0++)
		{
		  Boxes[count].resize(6);
		  int j0a=j0b+1;
		  j0b=((m0+1)*length)/FractalNodes0-1;
		  Boxes[count][0]=j0a;
		  Boxes[count][1]=j0b;
		  Boxes[count][2]=j1a;
		  Boxes[count][3]=j1b;
		  Boxes[count][4]=j2a;
		  Boxes[count][5]=j2b;
		  count++;
		}
	    }
	}
      assert(count==FractalNodes);
    }
    void calc_Buffers_and_more()
    {
      Buffers.resize(FractalNodes);
      BBoxes.resize(FractalNodes);
      PBoxes.resize(FractalNodes);
      PBoxesLength.resize(FractalNodes);
      BoxesLev.resize(FractalNodes);
      BBoxesLev.resize(FractalNodes);
      PBoxesLev.resize(FractalNodes);
      HRBoxesLev.resize(FractalNodes);
      HSBoxesLev.resize(FractalNodes);

      for(int count=0;count<FractalNodes;count++)
	{
	  Buffers[count].resize(6);
	  BBoxes[count].resize(6);
	  PBoxes[count].resize(6);
	  PBoxesLength[count].resize(3);
	  for(int n=0;n<3;n++)
	    {
	      if(Boxes[count][2*n] == 0 && !periodic)
		Buffers[count][2*n]=0;
	      else
		Buffers[count][2*n]=1;
	      if(Boxes[count][2*n+1] == grid_length && !periodic)
		Buffers[count][2*n+1]=0;
	      else
		Buffers[count][2*n+1]=1;
	      BBoxes[count][2*n]=Boxes[count][2*n]-Buffers[count][2*n];
	      BBoxes[count][2*n+1]=Boxes[count][2*n+1]+Buffers[count][2*n+1];
	      PBoxes[count][2*n]=BBoxes[count][2*n]-Buffers[count][2*n];
	      PBoxes[count][2*n+1]=BBoxes[count][2*n+1]+Buffers[count][2*n+1];
	      PBoxesLength[count][n]=PBoxes[count][2*n+1]-PBoxes[count][2*n]+1;
	    }
	  BoxesLev[count].resize(level_max+1);
	  BBoxesLev[count].resize(level_max+1);
	  PBoxesLev[count].resize(level_max+1);
	  HRBoxesLev[count].resize(level_max+1);
	  HSBoxesLev[count].resize(level_max+1);
	  BoxesLev[count][0].resize(6);
	  BBoxesLev[count][0].resize(6);
	  PBoxesLev[count][0].resize(6);
	  int zoom=Misc::pow(2,level_max);
	  for(int n=0;n<6;n++)
	    {
	      BoxesLev[count][0][n]=Boxes[count][n]*zoom;
	      BBoxesLev[count][0][n]=BBoxes[count][n]*zoom;
	      PBoxesLev[count][0][n]=PBoxes[count][n]*zoom;
	    }
	  for(int lev=1;lev<=level_max;lev++)
	    {
	      BoxesLev[count][lev].resize(6);
	      BBoxesLev[count][lev].resize(6);
	      PBoxesLev[count][lev].resize(6);
	      PBoxesLev[count][lev].resize(6);
	      HRBoxesLev[count][lev].resize(6);
	      HSBoxesLev[count][lev].resize(6);
	      zoom=Misc::pow(2,level_max-lev);
	      for(int n=0;n<3;n++)
		{
		  BoxesLev[count][lev][2*n]=BoxesLev[count][lev-1][2*n];
		  BBoxesLev[count][lev][2*n+1]=BBoxesLev[count][lev-1][2*n+1];
		  
		  BoxesLev[count][lev][2*n+1]=BBoxesLev[count][lev][2*n+1]-zoom*Buffers[count][2*n+1];
		  BBoxesLev[count][lev][2*n]=BoxesLev[count][lev][2*n]-zoom*Buffers[count][2*n];
		  
		  PBoxesLev[count][lev][2*n+1]=BBoxesLev[count][lev][2*n+1]+zoom*Buffers[count][2*n+1];
		  PBoxesLev[count][lev][2*n]=BBoxesLev[count][lev][2*n]-zoom*Buffers[count][2*n];
		  
		  HRBoxesLev[count][lev][2*n+1]=PBoxesLev[count][lev][2*n+1];
		  HRBoxesLev[count][lev][2*n]=BBoxesLev[count][lev][2*n];
		  
		  HSBoxesLev[count][lev][2*n+1]=BoxesLev[count][lev][2*n+1];
		  HSBoxesLev[count][lev][2*n]=BoxesLev[count][lev][2*n]+zoom*Buffers[count][2*n];
		}
	    }
	}
    }
    void calc_RealBoxes()
    {
      cout << "real " << FractalNodes << " " << grid_length << endl;
      RealBoxes.resize(FractalNodes);
      RealPBoxes.resize(FractalNodes);
      double glinv=1.0/static_cast<double>(grid_length);
      for(int b=0;b<FractalNodes;b++)
	{
	  RealBoxes[b].resize(6);
	  RealPBoxes[b].resize(6);
	  for(int ni=0;ni<6;ni+=2)
	    {
	      cout << " b ni " << b << " " << ni << endl;
	      RealBoxes[b][ni]=static_cast<double>(Boxes[b][ni])*glinv;
	      RealBoxes[b][ni+1]=static_cast<double>(Boxes[b][ni+1]+1)*glinv;
	      RealPBoxes[b][ni]=static_cast<double>(PBoxes[b][ni])*glinv;
	      RealPBoxes[b][ni+1]=static_cast<double>(PBoxes[b][ni+1])*glinv;
	      if(periodic)
		continue;
	      RealBoxes[b][ni]=max(RealBoxes[b][ni],glinv);
	      RealBoxes[b][ni+1]=min(RealBoxes[b][ni+1],1.0-glinv);
	      RealPBoxes[b][ni]=max(RealPBoxes[b][ni],glinv);
	      RealPBoxes[b][ni+1]=min(RealPBoxes[b][ni+1],1.0-glinv);
	    }
	}
    cout << " real b " << endl;
    }
    int fftw_where(const int& i,const int& j,const int& k,const int& lb,const int& lc)
    {
      return k+(j+(i-p_mess->start_x)*lb)*lc;
    }
    void make_scaling()
    {
      scaling=1.0;
      if(spectrum_number==1)
	{
	  double a1=pow(46.9*omega_0*h*h,0.67)*(1.0+pow(32.1*omega_0*h*h,-0.532));
	  double a2=pow(12.0*omega_0*h*h,0.424)*(1.0+pow(45.0*omega_0*h*h,-0.582));
	  double alpha=pow(a1,-omega_b/omega_0)*pow(a2,-pow(omega_b/omega_0,3));
	  scaling=box_length*omega_0*h*h*sqrt(alpha)*pow(1.0-omega_b/omega_0,0.6);
	  cout << "scaling " << a1 << " " << a2 << " " << alpha << " " << " " << box_length << " " << h << " " << scaling << endl;
	}
    }
    // public interface functions
    void fractal_memory_setup(Fractal_Memory* PFM);
    void setBalance(int B);
    void addParticles(int first,int total,
		      vector <double>& xmin,vector <double>& xmax,
		      vector <double>& xpos,vector <double>& ypos,
		      vector <double>& zpos,vector <double>& masses);
    void getField(int first,int last,double G,
		  vector <double>& xmin,vector <double>& xmax,
		  vector <double>& pot,vector <double>& fx,
		  vector <double>& fy,vector <double>& fz);
    void setNumberParticles(int NP);
    void setFractalNodes(int FR0,int FR1,int FR2);
    void setFFTNodes(int FN);
    void setPeriodic(bool per);
    void setDebug(bool Db);
    void setGridLength(int GL);
    void setPadding(int PA);
    void setLevelMax(int LM);
    void setMinimumNumber(int MN);
    void setHypreIterations(int MHI);
    void setHypreTolerance(double HT);
    void setBaseDirectory(string BD);
    void setRunIdentifier(string RI);
    void setTimeTrial(bool tt);
    void setTalkToMe(MPI_Comm& ttm);
    // static functions
    static double hubble(const double& arad,const double& omega_0,const double& omega_lambda)
    {
      return sqrt(omega_0/pow(arad,3)+(1.0-omega_0-omega_lambda)/pow(arad,2)+omega_lambda);
    }
    static double omega (const double& arad,const double& omega_0, const double& omega_lambda)
    //
    {
      return omega_0/(pow(arad,3)*pow(hubble(arad,omega_0,omega_lambda),2));  
    }
    //
    static double lambda(const double& arad,const double& omega_0, const double& omega_lambda)
    //
    {
      return omega_lambda/pow(hubble(arad,omega_0,omega_lambda),2);  
    }
  };
}
#endif
