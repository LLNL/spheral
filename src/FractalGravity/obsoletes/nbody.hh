#ifndef _Nbody_Defined_
#define _Nbody_Defined_
namespace FractalSpace
{
  using namespace std;
  class Nbody
  {
  public:
    int steps;
    int number_steps_total;
    int number_steps_out;
    int number_particles;
    int grid_length;
    int moat_0;
    bool periodic;
    bool random_initial;
    int minimum_number;
    int level_max;
    int tweaks;
    int padding;
    double off_x;
    double off_y;
    double off_z;
    double force_max;
    double omega_0;
    double omega_b;
    double omega_lambda;
    double omega_start;
    double lambda_start;
    double h;
    double box_length;
    int norm_what;
    int spectrum_number;
    int highest_level_init;
    double power_slope;
    double norm_scale;
    double sigma_0;
    double step_length;
    double arad;
    double parad;
    double pexp;
    double time;
    double udda;
    double redshift_start;
    double potential_energy;
    double potential_energy_old;
    double kinetic_energy;
    double kinetic_energy_old;
    vector <double> pos_x;
    vector <double> pos_y;
    vector <double> pos_z;
    vector <double> vel_x;
    vector <double> vel_y;
    vector <double> vel_z;
    vector <double> particle_mass;
    //  vector <int> highest_level;
    int splits;
    vector <double> splits_center_x;
    vector <double> splits_center_y;
    vector <double> splits_center_z;
    vector <double> splits_rad_x;
    vector <double> splits_rad_y;
    vector <double> splits_rad_z;
    vector <bool> splits_square;
    vector <int> splits_level;
    int masks;
    vector <double> masks_center_x;
    vector <double> masks_center_y;
    vector <double> masks_center_z;
    vector <double> masks_rad_x;
    vector <double> masks_rad_y;
    vector <double> masks_rad_z;
    vector <bool> masks_square;
    vector <int> masks_level;
    Nbody()
    {
      // default values
      // replace with your own values
      steps=0;
      number_steps_total=503;
      number_steps_out=100;
      number_particles=262144;
      grid_length=64;
      moat_0=1;
      periodic=true;
      random_initial=false;
      minimum_number=8;
      level_max=8;
      tweaks=2;
      padding=0;
      off_x=0.5;
      off_y=0.5;
      off_z=0.5;
      force_max=-1.0;
      omega_0=1.0;
      omega_lambda=0.0;
      omega_start=-1.0;
      lambda_start=-1.0;
      h=0.7;
      norm_what=0;
      spectrum_number=0;
      highest_level_init=4;
      power_slope=-1.0;
      norm_scale=0.2;
      sigma_0=1.0;
      box_length=100.0;
      step_length=0.04;
      arad=1.0;
      pexp=0.77;
      time=2.0/3.0;
      potential_energy=0.0;
      potential_energy_old=0.0;
      kinetic_energy=0.0;
      kinetic_energy_old=0.0;
      udda=0.0;
      redshift_start=49;
      splits=2;
      splits_center_x.assign(splits,0.5);
      splits_center_y.assign(splits,0.5);
      splits_center_z.assign(splits,0.5);
      splits_rad_x.assign(splits,0.5);
      splits_rad_y.assign(splits,0.5);
      splits_rad_z.assign(splits,0.5);
      splits_level.assign(splits,0);
      splits_square.assign(splits,true);
      splits_center_x[0]=0.05;
      splits_center_y[0]=0.05;
      splits_center_z[0]=0.05;
      splits_rad_x[0]=0.1;
      splits_rad_y[0]=0.1;
      splits_rad_z[0]=0.1;
      splits_level[0]=1;
      splits_square[0]=false;
      splits_center_x[1]=0.05;
      splits_center_y[1]=0.05;
      splits_center_z[1]=0.05;
      splits_rad_x[1]=0.07;
      splits_rad_y[1]=0.07;
      splits_rad_z[1]=0.07;
      splits_level[1]=2;
      splits_square[1]=false;
      masks=1;
      masks_center_x.assign(masks,0.5);
      masks_center_y.assign(masks,0.5);
      masks_center_z.assign(masks,0.5);
      masks_rad_x.assign(masks,0.5);
      masks_rad_y.assign(masks,0.5);
      masks_rad_z.assign(masks,0.5);
      masks_level.assign(masks,0);
      masks_square.assign(masks,true);
      masks_center_x[0]=0.05;
      masks_center_y[0]=0.05;
      masks_center_z[0]=0.05;
      masks_rad_x[0]=0.25;
      masks_rad_y[0]=0.25;
      masks_rad_z[0]=0.25;
      masks_level[0]=4;
      masks_square[0]=true;
      assign_vectors(2,0.0);
    }
    ~Nbody()
    {}
    void assign_vectors(const int& n,const double& value)
    {
      pos_x.assign(n*number_particles,value);
      pos_y.assign(n*number_particles,value);
      pos_z.assign(n*number_particles,value);
      particle_mass.assign(n*number_particles,value);
      //  highest_level.assign(n*number_particles,0);
      cout << "create " << pos_x.size() << endl;
    }
    void resize_vectors(const int& n)
    {
      pos_x.resize(n);
      pos_y.resize(n);
      pos_z.resize(n);
      vel_x.resize(n,0.0);
      vel_y.resize(n,0.0);
      vel_z.resize(n,0.0);
      particle_mass.resize(n);
      //  highest_level.resize(n,0);
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
