#ifndef _Particle_Defined_
#define _Particle_Defined_
namespace FractalSpace
{
  class Particle
  {
    vector <double> phase_space;
    // 0 1 2 x,y,z
    // 3 4 5 vx,vy,vz
    vector <double> field;
    // 0 1 2 3 pot,fx,fy,fz
    // 4 density
    // 5 rad_max
    double mass;
    int highest_level;
    Group* p_highest_level_group;
    int node_world;
    int number_world;
    bool real_particle;
  public:
    static int number_particles;
    Particle():
      highest_level(0),
      p_highest_level_group(0),
      node_world(-1),
      number_world(-1),
      real_particle(true)
    {
      phase_space.assign(3,0.0);
      field.assign(4,0.0);
      number_particles++;
    }
    Particle(const int& i,const int& j):
      highest_level(0),
      p_highest_level_group(0),
      node_world(-1),
      number_world(-1),
      real_particle(true)
    {
      phase_space.resize(i);
      field.resize(j);
      number_particles++;
    }
    Particle(const int& i,const double& x,const int& j,const double& y):
      highest_level(0),
      p_highest_level_group(0),
      node_world(-1),
      number_world(-1),
      real_particle(true)
    {
      phase_space.assign(i,x);
      field.assign(j,y);
      number_particles++;
    }
    Particle(Particle& p,vector <double>& shift)
    {
      phase_space=p.phase_space;
      phase_space[0]+=shift[0];
      phase_space[1]+=shift[1];
      phase_space[2]+=shift[2];
      field=p.field;
      mass=p.mass;
      highest_level=p.highest_level;
      p_highest_level_group=p.p_highest_level_group;
      node_world=p.node_world;
      number_world=p.number_world;
      real_particle=p.real_particle;
    }
    ~Particle()
    {
      assert(this);
      number_particles--;
    }
    void set_world(const int& node_w,const int& num_w);
    void get_world(int& node_w,int& num_w) const;
    double get_mass() const;
    int get_highest_level() const;
    void set_highest_level(const int& lev);
    double get_r2(const double& x,const double& y,const double& z) const;
    double get_r(const double& x,const double& y,const double& z) const;
    void get_phase_field_sizes(int& nphase,int& nfield) const;
    void get_pos(vector <double>& pos) const;
    void get_posm(vector <double>& pos,double& m) const;
    void get_pos(double& posx,double& posy,double& posz) const;
    void get_vel(vector <double>& vel) const;
    void get_force(vector <double>& force) const;
    void get_phase(vector <double>& pos,vector <double>& vel) const;
    void get_field(vector <double>& pos,vector <double>& vel,vector <double>& force) const;
    void get_field_pf(vector <double>& pf) const;
    double get_potential() const;
    double get_density() const;
    double get_rad_max() const;
    Group* get_p_highest_level_group() const;
    void set_mass(const double& m);
    void set_pos(vector <double>& pos);
    void set_posm(vector <double>& pos,double& m);
    void set_posmIFR(const double& posx,const double& posy,const double& posz,const double& m,const int& num_w,const int& node_w);
    void set_vel(vector <double>& vel);
    void set_pos(const double& posx,const double& posy,const double& posz);
    void set_phase(vector <double>& pos,vector <double>& vel);
    void set_field(vector <double>& pos,vector <double>& vel,vector <double>& force);
    void set_field_pf(vector <double>& sum_pf);
    void add_field_pf(vector <double>& sum_pf);
    void set_field_pf(vector <double>& sum_pf,const double& scale);
    void set_field_pf(const double& p,const double& f1,const double& f2,const double& f3);
    void set_field_pf(const double& x);
    void set_density(const double& d);
    void set_potential(const double& pot);
    void set_rad_max(const double& r);
    void set_p_highest_level_group(Group* g);
    void field_resize(const int& j);
    void space_resize(const int& i);
    void set_real_particle(const bool& rp);
    bool get_real_particle() const;
    void dump(ofstream& FILE) const;
    template <typename T> void dump(ofstream& FILE,vector <T>& pott,vector <T>& f_x,vector <T>& f_y,vector <T>& f_z) const;
  };
}
#endif
