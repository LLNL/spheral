#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
//
namespace FractalSpace
{
  void Particle::set_world(const int& node_w,const int& num_w)
  {
    node_world=node_w;
    number_world=num_w;
  }
  void Particle::get_world(int& node_w,int& num_w) const
  {
    node_w=node_world;
    num_w=number_world;
  }
  double Particle::get_mass() const
  {
    return mass;
  }
  int Particle::get_highest_level() const
  {
    return highest_level;
  }
  void Particle::set_highest_level(const int& lev)
  {
    highest_level=lev;
  }
  double Particle::get_r2(const double& x,const double& y,const double& z) const
  {
    return pow(phase_space[0]-x,2)+pow(phase_space[1]-y,2)+pow(phase_space[2]-z,2);
  }
  double Particle::get_r(const double& x,const double& y,const double& z) const
  {
    return sqrt(get_r2(x,y,z));
  }    
  void Particle::get_phase_field_sizes(int& nphase,int& nfield) const
  {
    nphase=phase_space.size();
    nfield=field.size();
  }
  void Particle::get_pos(vector <double>& pos) const
  {
    assert(phase_space.size() >= 3);
    pos[0]=phase_space[0];
    pos[1]=phase_space[1];
    pos[2]=phase_space[2];
  }
  void Particle::get_posm(vector <double>& pos,double& m) const
  {
    assert(phase_space.size() >= 3);
    pos[0]=phase_space[0];
    pos[1]=phase_space[1];
    pos[2]=phase_space[2];
    m=mass;
  }
  void Particle::get_pos(double& posx,double& posy,double& posz) const
  {
    assert(phase_space.size() >= 3);
    posx=phase_space[0];
    posy=phase_space[1];
    posz=phase_space[2];
  }
  void Particle::get_vel(vector <double>& vel) const
  {
    assert(phase_space.size() >= 6);
    vel[0]=phase_space[3];
    vel[1]=phase_space[4];
    vel[2]=phase_space[5];
  }
  void Particle::get_force(vector <double>& force) const
  {
    assert(field.size() >= 4);
    force[0]=field[1];
    force[1]=field[2];
    force[2]=field[3];
  }
  void Particle::get_phase(vector <double>& pos,vector <double>& vel) const
  {
    assert(phase_space.size() >= 6);
    pos[0]=phase_space[0];
    pos[1]=phase_space[1];
    pos[2]=phase_space[2];
    vel[0]=phase_space[3];
    vel[1]=phase_space[4];
    vel[2]=phase_space[5];
  }
  void Particle::get_field(vector <double>& pos,vector <double>& vel,vector <double>& force) const
  {
    assert(phase_space.size() >= 6);
    assert(field.size() >= 4);
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
  void Particle::get_field_pf(vector <double>& pf) const
  {
    assert(field.size() >= 4);
    pf[0]=field[0];
    pf[1]=field[1];
    pf[2]=field[2];
    pf[3]=field[3];
  }
  double Particle::get_potential() const
  {
    return field[0];
  }
  double Particle::get_density() const
  {
    return field[4];
  }
  double Particle::get_rad_max() const
  {
    return field[5];
  }
  Group* Particle::get_p_highest_level_group() const
  {
    return p_highest_level_group;
  }
  void Particle::set_mass(const double& m)
  {
    mass=m;
  }
  void Particle::set_pos(vector <double>& pos)
  {
    assert(phase_space.size() >= 3);
    phase_space[0]=pos[0];
    phase_space[1]=pos[1];
    phase_space[2]=pos[2];
  }
  void Particle::set_posm(vector <double>& pos,double& m) 
  {
    assert(phase_space.size() >= 3);
    phase_space[0]=pos[0];
    phase_space[1]=pos[1];
    phase_space[2]=pos[2];
    mass=m;
  }
  void Particle::set_posmIFR(const double& posx,const double& posy,const double& posz,const double& m,const int& num_w,const int& node_w)
  {
    assert(phase_space.size() >= 3);
    phase_space[0]=posx;
    phase_space[1]=posy;
    phase_space[2]=posz;
    mass=m;
    number_world=num_w;
    node_world=node_w;
  }
  void Particle::set_vel(vector <double>& vel)
  {
    assert(phase_space.size() >= 6);
    phase_space[3]=vel[0];
    phase_space[4]=vel[1];
    phase_space[5]=vel[2];
  }
  void Particle::set_pos(const double& posx,const double& posy,const double& posz)
  {
    assert(phase_space.size() >= 3);
    phase_space[0]=posx;
    phase_space[1]=posy;
    phase_space[2]=posz;
  }
  void Particle::set_phase(vector <double>& pos,vector <double>& vel)
  {
    assert(phase_space.size() >= 6);
    phase_space[0]=pos[0];
    phase_space[1]=pos[1];
    phase_space[2]=pos[2];
    phase_space[3]=vel[0];
    phase_space[4]=vel[1];
    phase_space[5]=vel[2];
  }
  void Particle::set_field(vector <double>& pos,vector <double>& vel,vector <double>& force)
  {
    assert(phase_space.size() >= 6);
    phase_space[0]=pos[0];
    phase_space[1]=pos[1];
    phase_space[2]=pos[2];
    phase_space[3]=vel[0];
    phase_space[4]=vel[1];
    phase_space[5]=vel[2];
    assert(field.size() >= 4);
    field[1]=force[0];
    field[2]=force[1];
    field[3]=force[2];
  }
  void Particle::set_field_pf(vector <double>& sum_pf)
  {
    assert(field.size() >= 4);
    field[0]=sum_pf[0];
    field[1]=sum_pf[1];
    field[2]=sum_pf[2];
    field[3]=sum_pf[3];
  }
  void Particle::add_field_pf(vector <double>& sum_pf)
  {
    assert(field.size() >= 4);
    field[0]+=sum_pf[0];
    field[1]+=sum_pf[1];
    field[2]+=sum_pf[2];
    field[3]+=sum_pf[3];
  }
  void Particle::set_field_pf(vector <double>& sum_pf,const double& scale)
  {
    assert(field.size() >= 4);
    field[0]=sum_pf[0]*scale;
    field[1]=sum_pf[1]*scale;
    field[2]=sum_pf[2]*scale;
    field[3]=sum_pf[3]*scale;
  }
  void Particle::set_field_pf(const double& p,const double& f1,const double& f2,const double& f3)
  {
    assert(field.size() >= 4);
    field[0]=p;
    field[1]=f1;
    field[2]=f2;
    field[3]=f3;
  }
  void Particle::set_field_pf(const double& x)
  {
    assert(field.size() >= 4);
    field[0]=x;
    field[1]=x;
    field[2]=x;
    field[3]=x;
  }
  void Particle::set_density(const double& d)
  {
    assert(field.size() > 4);
    field[4]=d;
  }
  void Particle::set_potential(const double& pot)
  {
    field[0]=pot;
  }
  void Particle::set_rad_max(const double& r)
  {
    assert(field.size() > 4);
    field[5]=r;
  }
  void Particle::set_p_highest_level_group(Group* g)
  {
    p_highest_level_group=g;
  }
  void Particle::field_resize(const int& j)
  {
    field.resize(j);
  }
  void Particle::space_resize(const int& i)
  {
    phase_space.resize(i);
  }
  void Particle::set_real_particle(const bool& rp)
  {
    real_particle=rp;
  }
  bool Particle::get_real_particle() const
  {
    return real_particle;
  }
  void Particle::dump(ofstream& FILE) const
  {
    FILE << " particle dump a ";
    FILE << "\t" << mass;
    FILE << "\t" << highest_level;
    FILE << "\t" << p_highest_level_group;
    FILE << "\t" << node_world;
    FILE << "\t" << number_world;
    FILE << "\t" << real_particle;
    FILE << "\t" << scientific;
    FILE << "\t" << phase_space[0];
    FILE << "\t" << phase_space[1];
    FILE << "\t" << phase_space[2];
    if(phase_space.size() == 6)
      {
	FILE << "\t" << phase_space[3];
	FILE << "\t" << phase_space[4];
	FILE << "\t" << phase_space[5];	  
      }
    FILE << "\t" << field[0];
    FILE << "\t" << field[1];
    FILE << "\t" << field[2];
    FILE << "\t" << field[3];
    if(field.size() > 4)
      {
	FILE << "\t" << field[4];
	FILE << "\t" << field[5];
      }
    else
      FILE << "\t" << " watch ";
    FILE << "\n";
  }
  template <typename T> void Particle::dump(ofstream& FILE,vector <T>& pott,vector <T>& f_x,vector <T>& f_y,vector <T>& f_z) const
  {
    FILE << " particle dump b " << "\n";
    FILE << scientific;
    FILE << "\t" << phase_space[0];
    FILE << "\t" << phase_space[1];
    FILE << "\t" << phase_space[2];
    FILE << "\t" << field[0];
    FILE << "\t" << field[1];
    FILE << "\t" << field[2];
    FILE << "\t" << field[3];
    FILE << "\n";
    for(int ii=0;ii < 8;ii++)
      FILE << pott[ii] << "\t" << f_x[ii] << "\t" << f_y[ii] << "\t" << f_z[ii] << "\n";
  }
  template void Particle::dump(ofstream& FILE,vector <double>& pott,vector <double>& f_x,vector <double>& f_y,vector <double>& f_z) const;
}
