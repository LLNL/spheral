#ifndef _Group_Defined_
#define _Group_Defined_

#include <vector>
#include <deque>

namespace FractalSpace
{
  class Group{
    bool done_group;
    bool buffer_group;
    int level;
    int group_number;
    int number_high_groups;
    Group* p_mother_group;
    int particles_in_group;
    int points_in_group;
    int number_high_points;
    int number_high_pairs;
    double force_const;
    double rjac;
    bool set_dens;
    bool set_scaling;
    Group* p_generated_from_group;
    std::vector <double> force_sum;
    double mass_sum;
    std::vector <int> miny;
    std::vector <int> maxy;
  public:
    File* p_file;
    std::vector <Point*> p_list_really_high;
    std::deque <Point*>list_points;
    std::vector <Point*>list_high_points;
    std::vector <Point*>list_new_points;
    std::vector <Group*>list_high_groups;
    std::vector <int> head_number;
    std::vector <Point*> list_high;
    std::vector <int> list_pair_1;
    std::vector <int> list_pair_2;
    static int number_groups;
    Group()
    {
      assert(this);
      done_group=false;
      buffer_group=true;
      p_file=0;
      level=0;
      group_number=-1;
      number_high_groups=-1;
      number_high_points=-1;
      number_high_pairs=-1;
      p_mother_group=0;
      particles_in_group=0;
      points_in_group=0;
      p_generated_from_group=this;
      force_sum.assign(3,0.0);
      mass_sum=1.0e-30;
      miny={INT_MAX,INT_MAX,INT_MAX};
      maxy={INT_MIN,INT_MIN,INT_MIN};
      number_groups++;
      //    cerr << "Making Group " << this << " " << number_groups << "\n";
    }
    Group(Group& mother_group)
    {
      assert(this);
      done_group=false;
      p_file=mother_group.p_file;
      buffer_group=mother_group.buffer_group;
      group_number=-1;
      level=mother_group.level+1;
      number_high_groups=-1;
      number_high_points=-1;
      number_high_pairs=-1;
      p_mother_group=&mother_group;
      particles_in_group=0;
      points_in_group=0;
      p_generated_from_group=this;
      force_sum.assign(3,0.0);
      mass_sum=1.0e-30;
      miny={INT_MAX,INT_MAX,INT_MAX};
      maxy={INT_MIN,INT_MIN,INT_MIN};
      number_groups++;
      //    cerr << "Making Group from Mother" << this << " " << number_groups << "\n";
    }
    ~Group()
    {
      number_groups--;
      //    cerr << "Ending Group " << this << " " << number_groups << "\n";
    }
    void set_group_number(const int& gn);
    int get_group_number() const;
    void set_forcem(const std::vector <double>& fs,const double& ms);
    void get_forcem(std::vector <double>& fs,double& ms) const;
    void set_level(const int& lev);
    int get_level() const;
    void set_force_const(const double& g_c);
    void set_rjac(const double& rj);
    double get_force_const() const;
    double get_rjac() const;
    Group* get_mother_group() const;
    void set_mother_group(Group* p_group);
    void set_points_in_group(const int& p);
    int get_points_in_group() const;
    void set_particles_in_group(const int& p);
    int get_number_high_points() const;
    void set_number_high_points(const int& i);
    int get_number_high_pairs() const;
    void set_number_high_pairs(const int& i);
    int get_number_high_groups() const;
    void set_number_high_groups(const int& i);
    void set_p_generated_from_group(Group* p_group);
    Group* get_p_generated_from_group() const;
    Group* get_p_mother_group() const;
    bool get_done_group() const;
    void set_done_group(bool b);
    bool get_buffer_group() const;
    void set_buffer_group(const bool& b);
    bool get_set_dens() const;
    void set_set_dens(const bool& b);
    bool get_set_scaling() const;
    void set_set_scaling(const bool& b);
    void subtract_density(const double& d);
    void scale_pot_forces(const double& scaling);
    void get_force_variance(double& varx,double& vary,double& varz);
    void get_miny_maxy(std::vector<int>& miny,std::vector<int>& maxy);
    void set_miny_maxy();
  };
}
#endif
