#ifndef _Fractal_Defined_
#define _Fractal_Defined_
namespace FractalSpace
{
  class Fractal
  {
    int number_particles;
    int grid_length;
    int length_ratio;
    bool periodic;
    bool halo;
    bool halo_fixed;
    unsigned int minimum_number;
    int level_max;
    int padding;
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
    int tweaks;
    vector <long> time_1;
    vector <long> time_2;
    vector <long>delta_time;
    vector <long> total_time;
    vector <long>time_g;
    vector <long>delta_g;
    vector <long>total_g;
    vector <long>time_p;
    vector <long>delta_p;
    vector <long>total_p;
    bool debug;
    int highest_level_used;
    int memory_value;
    Fractal_Memory* p_generated_from;
    vector <string> time_string;
    int steps;
  public:
    vector <Particle*> particle_list;
    vector <Particle*> particle_list_sorted;
    double omega_fraction;
    vector <double>rad;
    vector <double>grow;
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
      grid_length(64),
      length_ratio(1),
      periodic(true),
      halo(false),
      halo_fixed(false),
      minimum_number(8),
      level_max(8),
      padding(0),
      masks(0),
      epsilon_sor(6.0e-5),
      force_max(-1.0),
      halo_scale(1.0),
      halo_density0(1.0),
      density_0(0.0),
      moat_0(1),
      random_offset(0),
      maxits(250),
      tweaks(2),
      debug(false),
      highest_level_used(0),
      memory_value(0),
      omega_fraction(2.0/3.0)
    {
      steps=0;
      p_generated_from=0;
      assert(grid_length % length_ratio ==0);
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
      time_string[24]="Halo Force";
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
      memory_value(0),
      omega_fraction(2.0/3.0)
    {
      steps=0;
      p_generated_from=&mem;
      number_particles=mem.number_particles;
      grid_length=mem.grid_length;
      length_ratio=mem.length_ratio;
      periodic=mem.periodic;
      halo=mem.halo;
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
      tweaks=mem.tweaks;
      debug=mem.debug;
      assert(grid_length % length_ratio ==0);
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
      time_string[24]="Halo Force";
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
    int get_length_ratio()
    {
      return length_ratio;
    }
    void set_length_ratio(const int& i)
    {
      length_ratio=i;
      assert(grid_length % length_ratio==0);
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
    int get_tweaks()
    {
      return tweaks;
    }
    void set_tweaks(const int& i)
    {
      tweaks=i;
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
    void set_halo(const bool& i)
    {
      halo=i;
    }
    bool get_halo()
    {
      return halo;
    }
    bool get_halo_fixed()
    {
      return halo_fixed;
    }
    void set_memory_value(const int& i)
    {
      memory_value=i;
    }
    void copy_particle_list_to_list_sorted()
    {
      particle_list_sorted=particle_list;
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
    static double my_rand(const double& rand_max)
    {
      return (double)(rand())/rand_max;
    }
    static double my_rand_not_zero(const double& rand_max)
    {
      return (double)(max(rand(),1))/rand_max;
    }
  };
}
#endif
