#ifndef _Fractal_Memory_Defined_
#define _Fractal_memory_Defined_
namespace FractalSpace
{
  class Fractal_Memory
  {
  public:
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
    bool halo;
    bool halo_fixed;
    bool remember_points;
    int total_points_counter;
    int total_points_used;
    int total_points_generated;
    int new_points_gen;
    int random_gen;
    int highest_level_init;
    int norm_what;
    int spectrum_number;
    int particles_0_1d;
    int number_particles;
    int grid_length;
    int moat_0;
    unsigned int minimum_number;
    int tweaks;
    int padding;
    int level_max;
    int length_ratio;
    int number_steps_total;
    int number_steps_out;
    int random_offset;
    int maxits;
    unsigned int hypre_minimum;
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
    int node_points_max;
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
    Fractal_Memory():
      // default values
      // replace with your own values
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
      halo(false),
      halo_fixed(false),
      remember_points(false),
      total_points_counter(0),
      total_points_used(0),
      total_points_generated(0),
      new_points_gen(9),
      random_gen(12345),
      highest_level_init(4),
      norm_what(0),
      spectrum_number(0),
      particles_0_1d(64),
      number_particles(0),
      grid_length(64),
      moat_0(1),
      minimum_number(8),
      tweaks(2),
      padding(0),
      level_max(8),
      length_ratio(1),
      number_steps_total(503),
      number_steps_out(100),
      random_offset(0),
      maxits(250),
      hypre_minimum(1000),
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
      time(0.66666666667),
      total_mass(1.0),
      steps(0),
      //
      node_points_max(10000000),
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
    }
    ~Fractal_Memory()
    {
      cout << "Ending Fractal_Memory " << this << endl;
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
