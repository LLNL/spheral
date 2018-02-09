#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
//
namespace FractalSpace
{
  //
  void Group::set_group_number(const int& gn)
  {
    group_number=gn;
  }
  int Group::get_group_number() const
  {
    return group_number;
  }
  void Group::set_forcem(const vector <double>& fs,const double& ms)
  {
    force_sum=fs;
    mass_sum=ms;
  }
  void Group::get_forcem(vector <double>& fs,double& ms) const
  {
    fs=force_sum;
    ms=mass_sum;
  }
  void Group::set_level(const int& lev)
  {
    level=lev;
  }
  int Group::get_level() const
  {
    return level;
  }
  void Group::set_force_const(const double& g_c)
  {
    force_const=g_c;
  }
  void Group::set_rjac(const double& rj)
  {
    rjac=rj;
  }
  double Group::get_force_const() const
  {
    return force_const;
  }
  double Group::get_rjac() const
  {
    return rjac;
  }
  Group* Group::get_mother_group() const
  {
    return p_mother_group;
  }
  void Group::set_mother_group(Group* p_group)
  {
    p_mother_group=p_group;
  }
  void Group::set_points_in_group(const int& p)
  {
    points_in_group=p;
  }
  int Group::get_points_in_group() const
  {
    return points_in_group;
  }
  void Group::set_particles_in_group(const int& p)
  {
    particles_in_group=p;
  }
  int Group::get_number_high_points() const
  {
    return number_high_points;
  }
  void Group::set_number_high_points(const int& i)
  {
    number_high_points=i;
  }
  int Group::get_number_high_pairs() const
  {
    return number_high_pairs;
  }
  void Group::set_number_high_pairs(const int& i)
  {
    number_high_pairs=i;
  }
  int Group::get_number_high_groups() const
  {
    return number_high_groups;
  }
  void Group::set_number_high_groups(const int& i)
  {
    number_high_groups=i;
  }
  void Group::set_p_generated_from_group(Group* p_group)
  {
    p_generated_from_group=p_group;
  }
  Group* Group::get_p_generated_from_group() const
  {
    return p_generated_from_group;
  }
  Group* Group::get_p_mother_group() const
  {
    return p_mother_group;
  }
  bool Group::get_done_group() const
  {
    return done_group;
  }
  void Group::set_done_group(bool b)
  {
    done_group=b;
  }
  bool Group::get_buffer_group() const
  {
    return buffer_group;
  }
  void Group::set_buffer_group(const bool& b)
  {
    buffer_group=b;
  }
  bool Group::get_set_dens() const
  {
    return set_dens;
  }
  void Group::set_set_dens(const bool& b)
  {
    set_dens=b;
  }
  bool Group::get_set_scaling() const
  {
    return set_scaling;
  }
  void Group::set_set_scaling(const bool& b)
  {
    set_scaling=b;
  }
  void Group::subtract_density(const double& d)
  {
    for(auto p : list_points)
      p->subtract_dens_at_point(d);
    // for(vector <Point*>::const_iterator point_itr=list_points.begin();point_itr != list_points.end();++point_itr)
      //	if((*point_itr)->get_inside()) 
      // (*point_itr)->subtract_dens_at_point(d);
  }
  void Group::scale_pot_forces(const double& scaling)
  {
    for(auto p : list_points)
      p->scale_pot_forces(scaling);
    // for(vector <Point*>::const_iterator point_itr=list_points.begin();point_itr != list_points.end();++point_itr)
    //   (*point_itr)->scale_pot_forces(scaling);
  }
  void Group::get_force_variance(double& varx,double& vary,double& varz)
  {
    varx=0.0;
    vary=0.0;
    varz=0.0;
    double sum0(1.0e-30);
    double sumx(0.0);
    double sumy(0.0);
    double sumz(0.0);
    vector <double> f(3);
    for(auto p : list_points)
    // for(vector <Point*>::const_iterator point_itr=list_points.begin();point_itr != list_points.end();++point_itr)
      {
	sum0+=1.0;
	p->get_force_point(f);
	// (*point_itr)->get_force_point(f);
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
  void Group::set_miny_maxy()
  {
    if(list_points.empty())
      {
	miny={INT_MAX,INT_MAX,INT_MAX};
	maxy={INT_MIN,INT_MIN,INT_MIN};
	return;
      }
    vector<int>pos(3);
    list_points[0]->get_pos_point(pos);
    maxy=pos;
    miny=maxy;
    for(auto p : list_points)
      {
	p->get_pos_point(pos);
	for(int ni : {0,1,2})
	  {
	    miny[ni]=min(miny[ni],pos[ni]);
	    maxy[ni]=max(maxy[ni],pos[ni]);
	  }
      }
  }
  void Group::get_miny_maxy(vector<int>& MINY,vector<int>& MAXY)
  {
    if(miny[0] > maxy[0])
      set_miny_maxy();
    MINY=miny;
    MAXY=maxy;
  }
}
