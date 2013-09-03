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
    bool buffer_particle;
  public:
    static int number_particles;
    Particle()
    {
      assert(this);
      phase_space.reserve(6);
      field.reserve(6);
      highest_level=0;
      p_highest_level_group=0;
      buffer_particle=false;
      number_particles++;
    }
    Particle(const bool& b_p)
    {
      assert(this);
      phase_space.reserve(6);
      field.reserve(6);
      highest_level=0;
      p_highest_level_group=0;
      buffer_particle=b_p;
      number_particles++;
    }
    Particle(const int& i,const int& j)
    {
      assert(this);
      phase_space.reserve(6);
      field.reserve(6);
      phase_space.resize(i);
      field.resize(j);
      highest_level=0;
      p_highest_level_group=0;
      buffer_particle=false;
      number_particles++;
    }
    Particle(const int& i,const int& j,const bool& b_p)
    {
      assert(this);
      phase_space.reserve(6);
      field.reserve(6);
      phase_space.resize(i);
      field.resize(j);
      highest_level=0;
      p_highest_level_group=0;
      buffer_particle=b_p;
      number_particles++;
    }
    Particle(const int& i,const double& x,const int& j,const double& y)
    {
      assert(this);
      phase_space.assign(i,x);
      field.assign(j,y);
      highest_level=0;
      p_highest_level_group=0;
      buffer_particle=false;
      number_particles++;
    }
    Particle(const int& i,const double& x,const int& j,const double& y,const bool& b_p)
    {
      assert(this);
      phase_space.assign(i,x);
      field.assign(j,y);
      highest_level=0;
      p_highest_level_group=0;
      buffer_particle=b_p;
      number_particles++;
    }
    ~Particle()
    {
      assert(this);
      number_particles--;
    }
    double get_mass()
    {
      return mass;
    }
    void set_buffer_particle(const bool& b_p)
    {
      buffer_particle=b_p;
    }
    bool get_buffer_particle()
    {
      return buffer_particle;
    }
    int get_highest_level()
    {
      return highest_level;
    }
    void set_highest_level(const int& lev)
    {
      highest_level=lev;
    }
    double get_r2(const double& x,const double& y,const double& z)
    {
      return pow(phase_space[0]-x,2)+pow(phase_space[1]-y,2)+pow(phase_space[2]-z,2);
    }
    double get_r(const double& x,const double& y,const double& z)
    {
      return sqrt(get_r2(x,y,z));
    }    
    void get_pos(vector <double>& pos)
    {
      pos[0]=phase_space[0];
      pos[1]=phase_space[1];
      pos[2]=phase_space[2];
    }
    void get_posm(vector <double>& pos,double& m)
    {
      pos[0]=phase_space[0];
      pos[1]=phase_space[1];
      pos[2]=phase_space[2];
      m=mass;
    }
    void get_pos(double& posx,double& posy,double& posz)
    {
      posx=phase_space[0];
      posy=phase_space[1];
      posz=phase_space[2];
    }
    void get_vel(vector <double>& vel)
    {
      vel[0]=phase_space[3];
      vel[1]=phase_space[4];
      vel[2]=phase_space[5];
    }
    void get_force(vector <double>& force)
    {
      force[0]=field[1];
      force[1]=field[2];
      force[2]=field[3];
    }
    void get_phase(vector <double>& pos,vector <double>& vel)
    {
      pos[0]=phase_space[0];
      pos[1]=phase_space[1];
      pos[2]=phase_space[2];
      vel[0]=phase_space[3];
      vel[1]=phase_space[4];
      vel[2]=phase_space[5];
    }
    void get_field(vector <double>& pos,vector <double>& vel,vector <double>& force)
    {
      pos[0]=phase_space[0];
      pos[1]=phase_space[1];
      pos[2]=phase_space[2];
      vel[0]=phase_space[3];
      vel[1]=phase_space[4];
      vel[2]=phase_space[5];
      force[0]=field[1];
      force[1]=field[2];
      force[2]=field[3];
    }
    void get_field_pf(vector <double>& sum_pf)
    {
      sum_pf[0]=field[0];
      sum_pf[1]=field[1];
      sum_pf[2]=field[2];
      sum_pf[3]=field[3];
    }
    double get_potential()
    {
      return field[0];
    }
    double get_density()
    {
      return field[4];
    }
    double get_rad_max()
    {
      return field[5];
    }
    Group* get_p_highest_level_group()
    {
      return p_highest_level_group;
    }
    void set_mass(const double& m)
    {
      mass=m;
    }
    void set_pos(vector <double>& pos)
    {
      phase_space[0]=pos[0];
      phase_space[1]=pos[1];
      phase_space[2]=pos[2];
    }
    void set_vel(vector <double>& vel)
    {
      phase_space[3]=vel[0];
      phase_space[4]=vel[1];
      phase_space[5]=vel[2];
    }
    void set_pos(const double& posx,const double& posy,const double& posz)
    {
      phase_space[0]=posx;
      phase_space[1]=posy;
      phase_space[2]=posz;
    }
    void set_phase(vector <double>& pos,vector <double>& vel)
    {
      phase_space[0]=pos[0];
      phase_space[1]=pos[1];
      phase_space[2]=pos[2];
      phase_space[3]=vel[0];
      phase_space[4]=vel[1];
      phase_space[5]=vel[2];
    }
    void set_field(vector <double>& pos,vector <double>& vel,vector <double>& force)
    {
      phase_space[0]=pos[0];
      phase_space[1]=pos[1];
      phase_space[2]=pos[2];
      phase_space[3]=vel[0];
      phase_space[4]=vel[1];
      phase_space[5]=vel[2];
      field[1]=force[0];
      field[2]=force[1];
      field[3]=force[2];
    }
    void set_field_pf(vector <double>& sum_pf)
    {
      field[0]=sum_pf[0];
      field[1]=sum_pf[1];
      field[2]=sum_pf[2];
      field[3]=sum_pf[3];
    }
    void add_field_pf(vector <double>& sum_pf)
    {
      field[0]+=sum_pf[0];
      field[1]+=sum_pf[1];
      field[2]+=sum_pf[2];
      field[3]+=sum_pf[3];
    }
    void set_field_pf(vector <double>& sum_pf,const double& scale)
    {
      field[0]=sum_pf[0]*scale;
      field[1]=sum_pf[1]*scale;
      field[2]=sum_pf[2]*scale;
      field[3]=sum_pf[3]*scale;
    }
    void set_field_pf(const double& x)
    {
      field[0]=x;
      field[1]=x;
      field[2]=x;
      field[3]=x;
    }
    void set_density(const double& d)
    {
      field[4]=d;
    }
    void set_potential(const double& pot)
    {
      field[0]=pot;
    }
    void set_rad_max(const double& r)
    {
      field[5]=r;
    }
    void set_p_highest_level_group(Group* g)
    {
      p_highest_level_group=g;
    }
    void field_resize(const int& j)
    {
      field.resize(j);
    }
    void space_resize(const int& i)
    {
      phase_space.resize(i);
    }
    void dump()
    {
      cout << " particle dump a " << endl;
      cout << "\t" << phase_space[0];
      cout << "\t" << phase_space[1];
      cout << "\t" << phase_space[2];
      cout << "\t" << field[0];
      cout << "\t" << field[1];
      cout << "\t" << field[2];
      cout << "\t" << field[3];
      cout << endl;
    }
    template <typename T> void dump(vector <T>& pott,vector <T>& f_x,vector <T>& f_y,vector <T>& f_z)
    {
      cout << " particle dump b " << endl;
      cout << "\t" << phase_space[0];
      cout << "\t" << phase_space[1];
      cout << "\t" << phase_space[2];
      cout << "\t" << field[0];
      cout << "\t" << field[1];
      cout << "\t" << field[2];
      cout << "\t" << field[3];
      cout << endl;
      for(int ii=0;ii < 8;ii++)
	cout << pott[ii] << "\t" << f_x[ii] << "\t" << f_y[ii] << "\t" << f_z[ii] << endl;
    }
  };
}
#endif
