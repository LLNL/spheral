#ifndef _Fractal_Defined_
#define _Fractal_Defined_

#include <string>
#include <vector>
#include <deque>

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
    std::vector <int> Box;
    std::vector <int> BBox;
    std::vector <int> PBox;
    std::vector <double> RealBox;
    std::vector <int> PBoxLength; 
    std::vector < std::vector <int> >BoxLev;
    std::vector < std::vector <int> >BBoxLev;
    std::vector < std::vector <int> >PBoxLev;
    std::vector <int> Buffer;
    double clocks_per_sec;
    //
    int masks;
    std::vector <double> masks_center_x;
    std::vector <double> masks_center_y;
    std::vector <double> masks_center_z;
    std::vector <double> masks_rad_x;
    std::vector <double> masks_rad_y;
    std::vector <double> masks_rad_z;
    std::vector <bool> masks_square;
    std::vector <int> masks_level;
    double HTOL;
    double force_max;
    double halo_scale;
    double halo_density0;
    double density_0;
    int random_offset;
    int maxits;
    std::vector <double> time_1;
    std::vector <double> time_2;
    std::vector <double> delta_time;
    std::vector <double> total_time;
    std::vector <double> time_g;
    std::vector <double> delta_g;
    std::vector <double> total_g;
    std::vector <double> time_p;
    std::vector <double> delta_p;
    std::vector <double> total_p;
    bool debug;
    int highest_level_used;
    Fractal_Memory* p_generated_from;
    std::vector <std::string> time_string;
    int steps;
    double omega_start;
    double base_mass;
  public:
    Mess* p_mess;
    File* p_file;
    std::deque <Particle*> particle_list;
    std::deque <Particle*> particle_list_world;
    std::deque <Particle*> pseudo_particle_list;
    Particle* part_list_tmp;
    double omega_fraction;
    std::vector <double> rad;
    std::vector <double> grow;
    static bool first_time_solver;
    static std::string power_spec;
    static std::string integrator;
    static std::string energy_method;
    static std::string sim_parameters;
    static std::string force_fixed;
    static std::string particles;
    static std::string pot_solver;
    static std::string vel;
    Fractal()
    {
      //      cerr << " made ghost fractal" << "\n";
    }
    template <class M> Fractal(M& mem):
      density_0(0.0),
      debug(false),
      highest_level_used(0),
      omega_fraction(2.0/3.0)
    {
      //      cerr << " starting fractal " << "\n";
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
      HTOL=mem.HTOL;
      force_max=mem.force_max;
      halo_scale=mem.halo_scale;
      halo_density0=mem.halo_density0;
      random_offset=mem.random_offset;
      maxits=mem.maxits;
      debug=mem.debug;
      base_mass=mem.base_mass;
      //
      //      cerr << " fractal start a " << FractalNodes0 << " " << FractalNodes1 << " " << FractalNodes2 << " " << FractalNodes << "\n";
      p_mess=mem.p_mess;
      p_file=mem.p_file;
      MPIrun=mem.MPIrun;
      FractalNodes=mem.FractalNodes;
      FractalNodes0=mem.FractalNodes0;
      FractalNodes1=mem.FractalNodes1;
      FractalNodes2=mem.FractalNodes2;
      //      cerr << " fractal start b " << FractalNodes0 << " " << FractalNodes1 << " " << FractalNodes2 << " " << FractalNodes << "\n";
      //      cerr << " fractal start c " << p_mess << " " << p_file << "\n";
      FractalRank=get_FractalRank();
      //      cerr << FractalRank << "\n";
      assert(FractalRank<FractalNodes);
      Box=mem.Boxes[FractalRank];
      BBox=mem.BBoxes[FractalRank];
      PBox=mem.PBoxes[FractalRank];
      PBoxLength=mem.PBoxesLength[FractalRank];
      Buffer=mem.Buffers[FractalRank];
      //      BoxLev=mem.BoxesLev[FractalRank];
      //      BBoxLev=mem.BBoxesLev[FractalRank];
      //      PBoxLev=mem.PBoxesLev[FractalRank];
      BoxLev=mem.FRBoxesLev;
      BBoxLev=mem.FRBBoxesLev;
      PBoxLev=mem.FRPBoxesLev;
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
      //      cerr << "Making Fractal " << this << "\n";
    }
    ~Fractal()
    {    
      //      cerr << "Ending Fractal " << this << "\n";
    }
    void redo(Fractal_Memory* PFM);
    int get_FractalRank() const;
    int get_FractalNodes() const;
    void set_omega_start(int& omega);
    double get_omega_start() const;
    double get_base_mass() const;
    void set_base_mass(double bm);
    void setBox(std::vector <int>& B);
    void setBBox(std::vector <int>& BB);
    void getPBox(std::vector <int>& PB) const;
    void setPBox(std::vector <int>& PB);
    void setBuffer(std::vector <int>& BB);
    void getBox(std::vector <int>& B) const;
    void getBBox(std::vector <int>& BB) const;
    void getBBoxLev(std::vector <int>& BB,const int& level) const;
    void getPBoxLev(std::vector <int>& PB,const int& level) const;
    void getPBoxLength(std::vector <int>& PBL) const;
    void getBuffer(std::vector <int>& Bu) const;
    void getRealBox(std::vector <double>& RB) const;
    void assign_edge_buffer_passive(Point& point,const int& level,bool& edge,bool& buff,bool& pass,bool& really) const;
    void setMPIrun(const bool& Mr);
    bool getMPIrun() const;
    void inside_edge_buffer_pass(std::vector <int>& n,bool& inside,bool& edge,bool& buff,bool& pass) const;
    Fractal_Memory* get_p_generated_from() const;
    void set_p_generated_from (Fractal_Memory* p);
    double get_halo_scale() const;
    double get_halo_density0() const;
    int get_grid_length() const;
    void set_grid_length(const int& i);
    int get_highest_level_used() const;
    void set_highest_level_used(const int& i);
    double get_density_0() const;
    void set_density_0(const double& d);
    int get_number_particles() const;
    void set_number_particles(const int& i);
    int get_number_particles_world() const;
    void set_number_particles_world(const int& i);
    unsigned int get_minimum_number() const;
    void set_minimum_number(const int& i);
    int get_level_max() const;
    void set_level_max(const int& i);
    int get_padding() const;
    void set_padding(const int& i);
    void set_Hypre_TOL(const double& HT);
    double get_Hypre_TOL() const;
    int get_random_offset() const;
    void set_maxits(const int& m);
    int get_maxits() const;
    template <class M> void set_masks(M& mem);
    void set_pos_mask(const int& i,const double& x, const double& y, const double& z, const double& r1, const double& r2,const double& r3);
    double get_level_mask(const int& i) const;
    void set_level_mask(const int& i, const int& k);
    void get_mask(const int& m,std::vector <double>& cen,std::vector <double>& rad,bool& square) const;
    int get_number_masks() const;
    void set_number_masks(const int& i);
    bool get_debug() const;
    void set_force_max(const double& f_max);
    double get_force_max() const;
    void set_periodic(const bool& i);
    bool get_periodic() const;
    bool get_halo_fixed() const;
    void set_steps(int& s);
    void add_steps(int& s);
    int get_steps() const;
    void print_list(int what);
    void print_list_world(int what);
    double get_delta_time(int which) const;
    void get_total_times(std::vector <double>& TT) const;
    void timing_lev(const int& what,const int& level);
    void timing(const int& what, const int& which);
    void where_6(const int& i,const int& j,const int& k,std::vector <int>& Boxu) const;
    int where_1(const int& i,const int& j,const int& k) const;
    int where(const int& nx,const int& ny,const int& nz,std::vector <int>& boX,std::vector <int>& boXL) const;
    void wrap(const int& p);
    void wrap();
    void wrap(Particle* p);
    static double my_rand(const double& rand_max)
    {
      return (double)(rand())/rand_max;
    }
    static double my_rand_not_zero(const double& rand_max)
    {
      return (double)(std::max(rand(),1))/rand_max;
    }
    template <class M> static bool equal(std::vector <M>& a,std::vector <M>& b)
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
