#ifndef _Fractal_Memory_Defined_
#define _Fractal_Memory_Defined_

#include <string>
#include <vector>
#include <array>

namespace FractalSpace
{
  class Fractal_Memory
  {
  public:
    //
    std::string BaseDirectory;
    std::string RUN;
    MPI_Comm FractalWorld;
    bool standalone;
    const bool MPIrun;
    bool balance;
    int FractalNodes;
    int FractalNodes0;
    int FractalNodes1;
    int FractalNodes2;
    int FFTNodes;
    int SpecialNodes;
    bool time_trial;
    int min_hypre_group_size;
    bool hypre_load_balance;
    int hypre_min_node_load;
    int hypre_max_node_load;
    double hypre_multiplier;
    double G;
    std::vector <double>xmin;
    std::vector <double>xmax;
    std::vector <int> TouchWhichBoxes;
    std::vector <int> Touchy;
    std::vector < std::vector <int> > Boxes;
    std::vector < std::vector <int> > BBoxes;
    std::vector < std::vector <int> > PBoxes;
    std::vector < std::vector <int> > Buffers;
    std::vector < std::vector <int> > FRBoxesLev;
    std::vector < std::vector <int> > FRBBoxesLev;
    std::vector < std::vector <int> > FRPBoxesLev;
    std::vector < std::vector < std::vector <int> > > BoxesLev;
    std::vector < std::vector < std::vector <int> > > BBoxesLev;
    std::vector < std::vector < std::vector <int> > > PBoxesLev;
    std::vector < std::vector <int> > PBoxesLength;
    std::vector < std::vector <double> > RealBoxes;
    std::vector < std::vector <double> > RealPBoxes;
    std::vector < std::vector <double> > RealBBoxes;
    std::vector < std::vector <double> > RealIBoxes;
    std::vector <std::vector <int> > LeftCorners;
    std::vector <double> BigBox;
    std::string hypre_solver;
    std::string hypre_precond;
    int global_level_max;
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
    std::array<int,16>minimum_number_var;
    int padding;
    int level_max;
    int number_steps_total;
    int number_steps_out;
    int random_offset;
    int maxits;
    double base_mass;
    double HTOL;
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
    int level;
    //
    int crash_levels;
    double crash_pow;
    double density_crash;
    int max_particles;
    int splits;
    std::vector <double> splits_center_x;
    std::vector <double> splits_center_y;
    std::vector <double> splits_center_z;
    std::vector <double> splits_rad_x;
    std::vector <double> splits_rad_y;
    std::vector <double> splits_rad_z;
    std::vector <bool> splits_square;
    std::vector <int> splits_many;
    int masks;
    std::vector <double> masks_center_x;
    std::vector <double> masks_center_y;
    std::vector <double> masks_center_z;
    std::vector <double> masks_rad_x;
    std::vector <double> masks_rad_y;
    std::vector <double> masks_rad_z;
    std::vector <bool> masks_square;
    std::vector <int> masks_level;
    int masks_init;
    std::vector <double> masks_center_x_init;
    std::vector <double> masks_center_y_init;
    std::vector <double> masks_center_z_init;
    std::vector <double> masks_rad_x_init;
    std::vector <double> masks_rad_y_init;
    std::vector <double> masks_rad_z_init;
    std::vector <bool> masks_square_init;
    std::vector <int> masks_level_init;
    //
    std::vector < std::vector<Group*> > all_groups;

    Misc* p_misc; 
    Fractal* p_fractal;
    Mess* p_mess;
    File* p_file;
    static MPI_Comm FRACTAL_UNIVERSE;
    Fractal_Memory():
      //
      // default values
      // replace with your own values
      //
      BaseDirectory("FFRRAACCTTAALL/"),
      RUN("abc"),
      FractalWorld(MPI_COMM_WORLD),
      standalone(true),
      MPIrun(true),
      balance(1),
      FractalNodes(1),
      FractalNodes0(1),
      FractalNodes1(1),
      FractalNodes2(1),
      FFTNodes(1234567),
      SpecialNodes(0),
      time_trial(true),
      min_hypre_group_size(45),
      hypre_load_balance(false),
      hypre_min_node_load(45),
      hypre_max_node_load(40000),
      hypre_multiplier(2.0),
      G(1.0),
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
      HTOL(1.0e-7),
      force_max(-1.0),
      halo_scale(1.0),
      halo_density0(1.0),
      off_x(0.501),
      off_y(0.502),
      off_z(0.503),
      sigma_initial(-1.0),
      scaling(1.0),
      norm_scale(0.2),
      power_slope(-1.0),
      cut_off(1.0e6),
      cut_off_init(16.0),
      pexp(0.77),
      step_length(0.033),
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
      level(-10),
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
      minimum_number_var={{8,8,7,7,6,6,6,5,5,5,4,4,4,4,4,4}};
      hypre_solver="PCG";
      hypre_precond="PFMG";
      global_level_max=level_max;
      padding=std::min(padding,1);
      split_particles= force_max > 0.0;
      xmin.assign(3,0.0);
      xmax.assign(3,1.0);
    //
    }
    ~Fractal_Memory()
    {
      // cerr << "Ending Fractal_Memory " << this << "\n";
    }
    //
    void set_G(double Cavendish);
    void get_G(double& Cavendish) const;
    void set_xmin(std::vector <double>& xm);
    void set_xmax(std::vector <double>& xm);
    void get_xmin(std::vector <double>& xm) const;
    void get_xmax(std::vector <double>& xm) const;
    void calc_FractalNodes();
    void calc_Buffers_and_more();
    void calc_RealBoxes();
    int fftw_where(int i,int j,int k,int lb,int lc) const;
    void Full_Stop() const;
    void make_scaling();
    // public interface functions
    void fractal_memory_setup(Fractal_Memory* PFM);
    void setBalance(int B);
    void addParticles(int first,int total,
				      std::vector <double>& xmin,std::vector <double>& xmax,
				      std::vector <double>& xpos,std::vector <double>& ypos,
				      std::vector <double>& zpos,std::vector <double>& masses);
    void getField(int first,int last,double G,
				  std::vector <double>& xmin,std::vector <double>& xmax,
				  std::vector <double>& pot,std::vector <double>& fx,
				  std::vector <double>& fy,std::vector <double>& fz) const;
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
    void setBaseDirectory(std::string BD);
    void setRunIdentifier(std::string RI);
    void setTimeTrial(bool tt);
    void setTalkToMe(MPI_Comm& ttm);
  // static functions
  static double hubble(const double& arad,const double& omega_0,const double& omega_lambda)
  {
    return sqrt(omega_0/pow(arad,3)+(1.0-omega_0-omega_lambda)/pow(arad,2)+omega_lambda);
  }
  static double omega(const double& arad,const double& omega_0, const double& omega_lambda)
  {
    return omega_0/(pow(arad,3)*pow(hubble(arad,omega_0,omega_lambda),2));  
  }
  static double lambda(const double& arad,const double& omega_0, const double& omega_lambda)
  {
    return omega_lambda/pow(hubble(arad,omega_0,omega_lambda),2);  
  }
  };
}
#endif
