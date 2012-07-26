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
  public:
    static int number_particles;
    Particle():
      highest_level(0),
      p_highest_level_group(0),
      node_world(-1),
      number_world(-1)
    {
      phase_space.assign(3,0.0);
      number_particles++;
    }
    Particle(const int& i,const int& j):
      highest_level(0),
      p_highest_level_group(0),
      node_world(-1),
      number_world(-1)
    {
      phase_space.resize(i);
      field.resize(j);
      number_particles++;
    }
    Particle(const int& i,const double& x,const int& j,const double& y):
      highest_level(0),
      p_highest_level_group(0),
      node_world(-1),
      number_world(-1)
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
    }
    ~Particle()
    {
      assert(this);
      number_particles--;
    }
    void set_world(const int& node_w,const int& num_w)
    {
      node_world=node_w;
      number_world=num_w;
    }
    void get_world(int& node_w,int& num_w)
    {
      node_w=node_world;
      num_w=number_world;
    }
    double get_mass()
    {
      return mass;
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
    void get_field_pf(vector <double>& pf)
    {
      pf[0]=field[0];
      pf[1]=field[1];
      pf[2]=field[2];
      pf[3]=field[3];
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
    void set_posm(vector <double>& pos,double& m)
    {
      phase_space[0]=pos[0];
      phase_space[1]=pos[1];
      phase_space[2]=pos[2];
      mass=m;
    }
    void set_posmIFR(const double& posx,const double& posy,const double& posz,const double& m,const int& num_w,const int& node_w)
    {
      phase_space[0]=posx;
      phase_space[1]=posy;
      phase_space[2]=posz;
      mass=m;
      number_world=num_w;
      node_world=node_w;
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
    void set_field_pf(const double& p,const double& f1,const double& f2,const double& f3)
    {
      field[0]=p;
      field[1]=f1;
      field[2]=f2;
      field[3]=f3;
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
