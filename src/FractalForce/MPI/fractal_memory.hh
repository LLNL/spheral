#ifndef _Fractal_Memory_Defined_
#define _Fractal_Memory_Defined_
namespace FractalSpace
{
  class Fractal_Memory
  {
  public:
    //
    bool MPIrun;
    int FractalNodes;
    int FractalNodes0;
    int FractalNodes1;
    int FractalNodes2;
    vector < vector <int> > Boxes;
    vector < vector <int> > BBoxes;
    vector < vector <int> > PBoxes;
    vector < vector <int> > Buffers;
    vector < vector <int> > Box_to_Slices;
    vector < vector <int> > Box_from_Slices;
    vector < vector <int> > Slice_to_Boxes;
    vector < vector <int> > Slice_from_Boxes;
    vector < vector < vector <int> > > Box_to_Slices_Boxes;
    vector < vector < vector <int> > > Slice_to_Boxes_Boxes;
    vector < vector < vector <int> > > BoxesLev;
    vector < vector < vector <int> > > BBoxesLev;
    vector < vector < vector <int> > > PBoxesLev;
    vector < vector <int> > PBoxesLength;
    vector < vector <bool> > Periods;
    vector < vector <double> > RealBoxes;
    vector < vector <double> > RealPBoxes;
    Mess* p_mess;
    string hypre_solver;
    string hypre_precond;
    //
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
    //
    Fractal_Memory():
      //
      // default values
      // replace with your own values
      //
      MPIrun(false),
      FractalNodes(1),
      FractalNodes0(1),
      FractalNodes1(1),
      FractalNodes2(1),
      amnesia(true),
      mind_wipe(false),
      fixed_potential(false),
      calc_shear(false),
      start_up(true),
      calc_density_particle(true),
      do_vel(false),
      do_var(true), 
      periodic(true),
      random_initial(false),
      debug(false),
      halo_fixed(false),
      momentum_conserve(true),
      total_points_counter(0),
      total_points_used(0),
      total_points_generated(0),
      new_points_gen(9),
      random_gen(12345),
      highest_level_init(4),
      norm_what(0),
      spectrum_number(0),
      number_particles(262144),
      grid_length(64),
      moat_0(1),
      minimum_number(8),
      padding(0),
      level_max(8),
      number_steps_total(503),
      number_steps_out(100),
      random_offset(0),
      maxits(250),
      epsilon_sor(6.0e-5),
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
      step_length(0.04),
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
      hypre_solver="PCG";
      hypre_precond="SMG";
      //
      calc_FractalNodes();
      calc_Buffers_and_more();
      calc_RealBoxes();
      if(MPIrun)
	highest_level_init=0;     
    }
    ~Fractal_Memory()
    {
      cout << "Ending Fractal_Memory " << this << endl;
    }
    void calc_RealBoxes()
    {
      RealBoxes.resize(FractalNodes);
      RealPBoxes.resize(FractalNodes);
      double glinv=1.0/static_cast<double>(grid_length);
      for(int b=0;b<FractalNodes;b++)
	{
	  RealBoxes[b].resize(6);
	  RealPBoxes[b].resize(6);
	  for(int ni=0;ni<6;ni+=2)
	    {
	      RealBoxes[b][ni]=static_cast<double>(Boxes[b][ni])*glinv;
	      RealBoxes[b][ni+1]=static_cast<double>(Boxes[b][ni+1]+1)*glinv;
	      RealPBoxes[b][ni]=static_cast<double>(PBoxes[b][ni])*glinv;
	      RealPBoxes[b][ni+1]=static_cast<double>(PBoxes[b][ni+1])*glinv;
	    }
	}
    }
    void calc_FractalNodes()
    {
      FractalNodes=FractalNodes0*FractalNodes1*FractalNodes2;
      MPIrun=FractalNodes > 1;
      Boxes.resize(FractalNodes);
      int count=0;
      int j0b=-1;
      for(int m0=0;m0<FractalNodes0;m0++)
	{
	  int j0a=j0b+1;
	  j0b=((m0+1)*grid_length)/FractalNodes0-1;
	  int j1b=-1;
	  for(int m1=0;m1<FractalNodes1;m1++)
	    {
	      int j1a=j1b+1;
	      j1b=((m1+1)*grid_length)/FractalNodes1-1;
	      int j2b=-1;
	      for(int m2=0;m2<FractalNodes2;m2++)
		{
		  Boxes[count].resize(6);
		  int j2a=j2b+1;
		  j2b=((m2+1)*grid_length)/FractalNodes2-1;
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
      Periods.resize(FractalNodes);
      Buffers.resize(FractalNodes);
      BBoxes.resize(FractalNodes);
      PBoxes.resize(FractalNodes);
      PBoxesLength.resize(FractalNodes);
      BoxesLev.resize(FractalNodes);
      BBoxesLev.resize(FractalNodes);
      PBoxesLev.resize(FractalNodes);

      for(int count=0;count<FractalNodes;count++)
	{
	  Periods[count].resize(3);
	  Buffers[count].resize(6);
	  BBoxes[count].resize(6);
	  PBoxes[count].resize(6);
	  PBoxesLength[count].resize(3);

	  Periods[count][0]=periodic && Boxes[count][0] == 0 && Boxes[count][1] == grid_length-1;
	  Periods[count][1]=periodic && Boxes[count][2] == 0 && Boxes[count][3] == grid_length-1;
	  Periods[count][2]=periodic && Boxes[count][4] == 0 && Boxes[count][5] == grid_length-1;
	  for(int n=0;n<3;n++)
	    {
	      if(Periods[count][n] || (Boxes[count][2*n] == 0 && !periodic))
		Buffers[count][2*n]=0;
	      else
		Buffers[count][2*n]=1;
	      if(Periods[count][n] || (Boxes[count][2*n+1] == grid_length-1 && !periodic))
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
	      zoom=Misc::pow(2,level_max-lev);
	      for(int n=0;n<3;n++)
		{
		  BoxesLev[count][lev][2*n]=BoxesLev[count][lev-1][2*n];
		  BBoxesLev[count][lev][2*n+1]=BBoxesLev[count][lev-1][2*n+1];
		  
		  BoxesLev[count][lev][2*n+1]=BBoxesLev[count][lev][2*n+1]-zoom*Buffers[count][2*n+1];
		  BBoxesLev[count][lev][2*n]=BoxesLev[count][lev][2*n]-zoom*Buffers[count][2*n];
		  
		  PBoxesLev[count][lev][2*n+1]=BBoxesLev[count][lev][2*n+1]+zoom*Buffers[count][2*n+1];
		  PBoxesLev[count][lev][2*n]=BBoxesLev[count][lev][2*n]-zoom*Buffers[count][2*n];
		}
	    }
	}
    }
    void calc_Boxes_vs_Slices()
    {
      Box_to_Slices.resize(FractalNodes);
      Box_from_Slices.resize(FractalNodes);
      Slice_to_Boxes.resize(FractalNodes);
      Slice_from_Boxes.resize(FractalNodes);
      for(int b=0;b<FractalNodes;b++)
	{
	  for(int s=0;s<FractalNodes;s++)
	    {
	      if(p_mess->Slices[s][0] <= Boxes[b][1] && p_mess->Slices[s][1] >= Boxes[b][0])
		{
		  Box_to_Slices[b].push_back(s);
		  Slice_from_Boxes[s].push_back(b);
		}
	      if(p_mess->Slices[s][0] <= BBoxes[b][1] && p_mess->Slices[s][1] >= BBoxes[b][0])
		{
		  Slice_to_Boxes[s].push_back(b);
		  Box_from_Slices[b].push_back(s);
		}
	    }
	}
      Box_to_Slices_Boxes.resize(FractalNodes);
      Slice_to_Boxes_Boxes.resize(FractalNodes);
      for(int b=0;b<FractalNodes;b++)
	{
	  Box_to_Slices_Boxes[b].resize(FractalNodes);
	  int bs_size=Box_to_Slices[b].size();
	  for(int s=0;s<bs_size;s++)
	    {
	      int sl=Box_to_Slices[b][s];
	      Box_to_Slices_Boxes[b][sl].push_back(max(p_mess->Slices[sl][0],Boxes[b][0]));
	      Box_to_Slices_Boxes[b][sl].push_back(min(p_mess->Slices[sl][1],Boxes[b][1]));
	    }
	}
      for(int s=0;s<FractalNodes;s++)
	{
	  Slice_to_Boxes_Boxes[s].resize(FractalNodes);
	  int sb_size=Slice_to_Boxes[s].size();
	  for(int b=0;b<sb_size;b++)
	    {
	      int bs=Slice_to_Boxes[s][b];
	      Slice_to_Boxes_Boxes[s][bs].push_back(max(p_mess->Slices[s][0],BBoxes[bs][0]));
	      Slice_to_Boxes_Boxes[s][bs].push_back(min(p_mess->Slices[s][1],BBoxes[bs][1]));
	    }
	}
    }
    int fftw_where(const int& i,const int& j,const int& k,const int& la,const int& lb)
    {
      return k+(j+i*la)*lb;
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
