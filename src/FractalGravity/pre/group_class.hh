#ifndef _Group_Defined_
#define _Group_Defined_
namespace FractalSpace
{
  class Group{
    bool done_group;
    int level;
    int number_high_groups;
    Group* p_mother_group;
    int particles_in_group;
    int points_in_group;
    int number_high_points;
    int number_high_pairs;
    double grav_const;
    double rjac;
    bool set_dens;
    bool set_scaling;
    Group* p_generated_from_group;
  public:
    vector <Point*> p_list_really_high;
    vector <Point*>list_points;
    vector <Point*>list_high_points;
    vector <Point*>list_new_points;
    vector <Group*>list_high_groups;
    vector <int> head_number;
    vector <Point*> list_high;
    vector <int> list_pair_1;
    vector <int> list_pair_2;
    static int number_groups;
    Group(){
      assert(this);
      done_group=false;
      level=0;
      number_high_groups=-1;
      number_high_points=-1;
      number_high_pairs=-1;
      p_mother_group=0;
      particles_in_group=0;
      points_in_group=0;
      p_generated_from_group=this;
      number_groups++;
      //    cout << "Making Group " << this << " " << number_groups << endl;
    };
    Group(Group& mother_group){
      assert(this);
      done_group=false;
      level=mother_group.level+1;
      number_high_groups=-1;
      number_high_points=-1;
      number_high_pairs=-1;
      p_mother_group=&mother_group;
      particles_in_group=0;
      points_in_group=0;
      p_generated_from_group=this;
      number_groups++;
      //    cout << "Making Group from Mother" << this << " " << number_groups << endl;
    };
    ~Group()
    {
      number_groups--;
      //    cout << "Ending Group " << this << " " << number_groups << endl;
    };
    /*
      Group(const Group&){};
      Group& operator=(const Group&){};
    */
    //
    void set_level(const int& lev)
    {
      level=lev;
    }
    int get_level()
    {
      return level;
    }
    void set_grav_const(const double& g_c)
    {
      grav_const=g_c;
    }
    void set_rjac(double& rj)
    {
      rjac=rj;
    }
    double get_grav_const()
    {
      return grav_const;
    }
    double get_rjac()
    {
      return rjac;
    }
    Group* get_mother_group()
    {
      return p_mother_group;
    }
    void set_mother_group(Group* p_group)
    {
      p_mother_group=p_group;
    }
    void set_points_in_group(const int& p)
    {
      points_in_group=p;
    }
    int get_points_in_group()
    {
      return points_in_group;
    }
    void set_particles_in_group(const int& p)
    {
      particles_in_group=p;
    }
    int get_number_high_points()
    {
      return number_high_points;
    }
    void set_number_high_points(const int& i)
    {
      number_high_points=i;
    }
    int get_number_high_pairs()
    {
      return number_high_pairs;
    }
    void set_number_high_pairs(const int& i)
    {
      number_high_pairs=i;
    }
    int get_number_high_groups()
    {
      return number_high_groups;
    }
    void set_number_high_groups(const int& i)
    {
      number_high_groups=i;
    }
    void set_p_generated_from_group(Group* p_group)
    {
      p_generated_from_group=p_group;
    }
    Group* get_p_generated_from_group()
    {
      return p_generated_from_group;
    }
    Group* get_p_mother_group()
    {
      return p_mother_group;
    }
    bool get_done_group()
    {
      return done_group;
    }
    void set_done_group(const bool& b)
    {
      done_group=b;
    }
    bool get_set_dens()
    {
      return set_dens;
    }
    void set_set_dens(const bool& b)
    {
      set_dens=b;
    }
    bool get_set_scaling()
    {
      return set_scaling;
    }
    void set_set_scaling(const bool& b)
    {
      set_scaling=b;
    }
    void subtract_density(const double& d)
    {
      for(vector <Point*>:: const_iterator point_itr=list_points.begin();point_itr != list_points.end();++point_itr)
	if((*point_itr)->get_inside()) (*point_itr)->subtract_dens_at_point(d);
    }
    void scale_pot_forces(const double& scaling)
    {
      for(vector <Point*>::const_iterator point_itr=list_points.begin();point_itr != list_points.end();++point_itr)
	{
	  Point* p_point=*point_itr;
	  p_point->scale_pot_forces(scaling);
	}
    }
    void get_force_variance(double& varx,double& vary,double& varz)
    {
      varx=0.0;
      vary=0.0;
      varz=0.0;
      double sum0=1.0e-30;
      double sumx=0.0;
      double sumy=0.0;
      double sumz=0.0;
      vector <double> f(3);
      for(vector <Point*>::const_iterator point_itr=list_points.begin();point_itr != list_points.end();++point_itr)
	{
	  Point* p_point=*point_itr;
	  sum0+=1.0;
	  p_point->get_force_point(f);
	  sumx+=f[0];
	  varx+=f[0]*f[0];
	  sumy+=f[1];
	  vary+=f[1]*f[1];
	  sumz+=f[2];
	  varz+=f[2]*f[2];
	}
      varx=sqrt(varx/sum0-pow(sumx/sum0,2));
      vary=sqrt(vary/sum0-pow(sumy/sum0,2));
      varz=sqrt(varz/sum0-pow(sumz/sum0,2));
    }
  };
}
#endif
