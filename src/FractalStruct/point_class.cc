#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
//
namespace FractalSpace
{
  void Point::set_found_it(bool fi)
  {
    found_it=fi;
  }
  bool Point::get_found_it() const
  {
    return found_it;
  }
  void Point::set_mass_point(bool what)
  {
    mass_point=what;
  }
  void Point::set_mass_points(bool what)
  {
    mass_point=what;
    Point* p=get_point_ud_0(1);
    p->mass_point=what;
    p=p->get_point_ud_0(3);
    p->mass_point=what;
    p=p->get_point_ud_0(0);
    p->mass_point=what;
    p=p->get_point_ud_0(5);
    p->mass_point=what;
    p=p->get_point_ud_0(1);
    p->mass_point=what;
    p=p->get_point_ud_0(2);
    p->mass_point=what;
    p=p->get_point_ud_0(0);
    p->mass_point=what;
  }
  bool Point::get_mass_point() const
  {
    return mass_point;
  }
  bool Point::get_trouble() const
  {
    return in_trouble;
  }
  void Point::set_trouble(bool what)
  {
    in_trouble=what;
  }
  void Point::set_edge_buffer_passive_point(const bool& e,const bool& b,const bool& p)
  {
    edge_point=e;
    buffer_point=b;
    passive_point=p;
  }
  void Point::set_edge_buffer_passive_really_point(const bool& e,const bool& b,const bool& p,const bool& r)
  {
    edge_point=e;
    buffer_point=b;
    passive_point=p;
    really_passive=r;
  }
  bool Point::get_buffer_point() const
  {
    return buffer_point;
  }
  void Point::set_buffer_point(bool& b_p)
  {
    buffer_point=b_p;
  }
  bool Point::get_edge_point() const
  {
    return edge_point;
  }
  void Point::set_edge_point(bool& e_p)
  {
    edge_point=e_p;
  }
  int Point::get_number_in_list() const
  {
    return number_in_list;
  }
  void Point::set_number_in_list(const int& n)
  {
    number_in_list=n;
  }
  void Point::set_eureka_adj(const vector <bool>& eu)
  {
    eureka_adj=eu;
  }
  void Point::get_eureka_adj(vector <bool>& eu) const
  {
    eu=eureka_adj;
  }
  void Point::set_eureka_dau(const vector <bool>& eu)
  {
    eureka_dau=eu;
  }
  void Point::get_eureka_dau(vector <bool>& eu) const
  {
    eu=eureka_dau;
  }
  bool Point::get_eureka_dau(const int& ni) const
  {
    return eureka_dau[ni];
  }
  void Point::set_p_in_group(Group* p_g) 
  {
    p_in_group=p_g;
  }
  Group* Point::get_p_in_group() const
  {
    return p_in_group;
  }
  void Point::set_p_in_high_group(Group* p_g) 
  {
    p_in_high_group=p_g;
  }
  Group* Point::get_p_in_high_group()
  {
    return p_in_high_group;
  }
  void Point::get_pos_point(array <int,3>& pos) const
  {
    pos=pos_point;
  }
  void Point::get_pos_point(vector <int>& pos) const
  {
    pos.assign(pos_point.begin(),pos_point.begin()+3);
  }
  vector<int> Point::get_pos_point() const
  {
    vector<int>pp(3);
    pp.assign(pos_point.begin(),pos_point.begin()+3);
    return pp;
  }
  array<int,3> Point::get_pos_point_a() const
  {
    array<int,3>pp=pos_point;
    return pp;
  }
  void Point::set_pos_point(const array <int,3>& pos)
  {
    pos_point=pos;
  }
  void Point::set_pos_point(const vector <int>& pos)
  {
    for(int ni : {0,1,2})
      pos_point[ni]=pos[ni];
  }
  void Point::set_pos_point(const int& x,const int& y,const int& z)
  {
    pos_point[0]=x;
    pos_point[1]=y;
    pos_point[2]=z;
  }
  void Point::get_pos_point(int& x,int& y,int& z) const
  {
    x=pos_point[0];
    y=pos_point[1];
    z=pos_point[2];
  }
  int Point::get_pos_point_x() const
  {
    return pos_point[0];
  }
  void Point::set_pos_point_x(const int& i)
  {
    pos_point[0]=i;
  }
  int Point::get_pos_point_y() const
  {
    return pos_point[1];
  }
  void Point::set_pos_point_y(const int& i)
  {
    pos_point[1]=i;
  }
  int Point::get_pos_point_z() const
  {
    return pos_point[2];
  }
  void Point::set_pos_point_z(const int& i)
  {
    pos_point[2]=i;
  }
  int Point::get_pos_point(const int& i) const
  {
    return pos_point[i];
  }
  int Point::get_real_pointer() const
  {
    return real_pointer;
  }
  void Point::set_real_pointer(const int& i)
  {
    real_pointer=i;
  }
  // void Point::set_point_to_number(const int& i)
  // {
  //   point_to_number=i;
  // }
  // int Point::get_point_to_number() const
  // {
  //   return point_to_number;
  // }
  void Point::set_high_number(const int& i)
  {
    high_number=i;
  }
  int Point::get_high_number() const
  {
    return high_number;
  }
  bool Point::get_inside() const
  {
    return  inside;
  }
  void Point::set_inside(const bool& ins)
  {
    inside=ins;
  }
  bool Point::get_passive_point() const
  {
    return passive_point;
  }
  void Point::set_passive_point(const bool& value)
  {
    passive_point=value;
  }
  bool Point::get_really_passive() const
  {
    return really_passive;
  }
  void Point::set_really_passive(const bool& value)
  {
    really_passive=value;
  }
  // void Point::set_ij_number(const int& count)
  // {
  //   ij_number=count;
  // }
  // int Point::get_ij_number() const
  // {
  //   return ij_number;
  // }
  void Point::really_clear(vector <Point*>& die)
  {
    vector <Point*>temp;
    die.swap(temp);
  }
  // void Point::set_ij_neighbors()
  // {
  //   if(!inside || passive_point)
  //     {
  // 	ij_ud.clear();
  // 	return;
  //     }
  //   ij_ud.resize(6);
  //   for(int ni=0;ni<6;ni++)
  //     {
  // 	assert(!point_ud[ni]->really_passive);
  // 	ij_ud[ni]=point_ud[ni]->ij_number;
  //     }
  // }
  // void Point::set_ij_neighbors(vector <int>& Box)
  // {
  //   if(!inside)
  //     {
  // 	ij_ud.clear();
  // 	return;
  //     }
  //   if(pos_point[0] > Box[0] && pos_point[0] < Box[1] &&
  //      pos_point[1] > Box[2] && pos_point[1] < Box[3] &&
  //      pos_point[2] > Box[4] && pos_point[2] < Box[5])
  //     {
  // 	ij_ud.resize(6);
  // 	for(int ni=0;ni<6;ni++)
  // 	  ij_ud[ni]=get_point_ud_0(ni)->ij_number;
  // 	return;
  //     }
  //   ij_ud.clear();
  //   return;
  // }
  // void Point::get_ij_neighbors(vector <int>& ijud) const
  // {
  //   ijud=ij_ud;
  // }
  // int Point::get_ij_neighbors_size() const
  // {
  //   return ij_ud.size();
  // }
  // void Point::copy_ij_index(const int& ijc)
  // {
  //   ij_ud.resize(1);
  //   ij_ud[0]=ijc;
  // }
  // void Point::get_hypre_info(int& ij_index,vector <int>& ijud,double& rho,double& pot) const
  // {
  //   ij_index=ij_number;
  //   ijud=ij_ud;
  //   rho=density_point;
  //   pot=potential_point;
  // }
  bool Point::get_it_is_high() const
  {
    return it_is_high;
  }
  void Point::set_it_is_high(const bool& value)
  {
    it_is_high=value;
  }
  bool Point::get_it_is_really_high() const
  {
    return it_is_really_high;
  }
  void Point::set_it_is_really_high(const bool& value)
  {
    it_is_really_high=value;
  }
  void Point::set_passive_low()
  {
    it_is_high=it_is_high && !passive_point;
  }
  void Point::set_inside_high()
  {
    it_is_high=it_is_high && inside; 
  }
  double Point::get_potential_point() const
  {
    return potential_point;
  }
  void Point::set_potential_point(const double& d)
  {
    potential_point=d;
  }
  void Point::add_potential_point(const double& d)
  {
    potential_point+=d;
  }
  void Point::add_force_point_x(const double& d)
  {
    force_point[0]+=d;
  }
  void Point::add_force_point_y(const double& d)
  {
    force_point[1]+=d;
  }
  void Point::add_force_point_z(const double& d)
  {
    force_point[2]+=d;
  }
  double Point::get_density_point() const
  {
    return density_point;
  }
  void Point::get_force_point(vector <double>& force) const
  {
//     force=force_point;
    force.assign(force_point.begin(),force_point.begin()+3);
    // if(force.size() < 3)
    //   force.resize(3);
    // std::copy(force_point,force_point+3,force.begin());
  }
  void Point::set_force_point(vector <double>& force)
  {
//     force_point=force;
    for(int ni : {0,1,2})
      force_point[ni]=force[ni];
    // std::copy(force.begin(),force.begin()+3,force_point);
  }
  double Point::get_force_point_x() const
  {
    return force_point[0];
  }
  //
  double Point::get_force_point_y() const
  {
    return force_point[1];
  }
  //
  double Point::get_force_point_z() const
  {
    return force_point[2];
  }
  void Point::set_force_point_x(const double& f)
  {
    force_point[0]=f;
  }
  void Point::set_force_point_y(const double& f)
  {
    force_point[1]=f;
  }
  void Point::set_force_point_z(const double& f)
  {
    force_point[2]=f;
  }
//   void Point::set_FILE(ofstream* p_filE)
//   {
//     p_FILE=p_filE;
//   }
//   ofstream* Point::get_FILE()
//   {
//     return p_FILE;
//   }
//   void Point::all_mine(vector <Point*>& pointers,vector <bool>& belongs_to_me)
//   {
//     assert(real_pointer==0);
//     belongs_to_me.resize(27);
//     pointers.resize(27);
//     int n=0;
//     Point* p_point_z=this;
//     for(int nz=0;nz<3;++nz)
//       {
// 	Point* p_point_y=p_point_z;
// 	for(int ny=0;ny<3;++ny)
// 	  {
// 	    Point* p_point_x=p_point_y; 
// 	    for(int nx=0;nx<3;++nx)
// 	      {
// 		pointers[n]=p_point_x;
// 		belongs_to_me[n]=p_point_x->real_pointer==n;
// 		if(nx < 2)
// 		  p_point_x=p_point_x->get_point_up_x_0();
// 		n++;
// 	      }
// 	    if(ny < 2)
// 	      p_point_y=p_point_y->get_point_up_y_0();
// 	  }
// 	if(nz < 2)
// 	  p_point_z=p_point_z->get_point_up_z_0();
//       }
//   }
  Point* Point::move_adj(const int& ra,const int& rb)
  {
    Point* arthur=this;
    int dr=rb/4 - ra/4;
    if(dr > 0) 
      arthur=arthur->get_point_up_z_0();
    else if(dr < 0) 
      arthur=arthur->get_point_down_z_0();
    dr=(rb/2 % 2) - (ra/2 % 2);
    if(dr > 0) 
      arthur=arthur->get_point_up_y_0();
    else if(dr < 0) 
      arthur=arthur->get_point_down_y_0();
    dr=(rb % 2) - (ra % 2);
    if(dr > 0) 
      arthur=arthur->get_point_up_x_0();
    else if(dr < 0) 
      arthur=arthur->get_point_down_x_0();
    return arthur;
  }
  Point* Point::move_rp(const int& r)
  {
    Point* arthur=this;
    int drx=(r % 3)-(real_pointer % 3);
    int dry=(r/3 %3)-(real_pointer/3 % 3);
    int drz=(r/9)-(real_pointer/9);
    if(drx > 0)
      {
	for (int dr=0;dr < drx;dr++)
	  arthur=arthur->get_point_ud_0(1);
      }
    else if(drx < 0)
      {
	for (int dr=0;dr < -drx;dr++)
	  arthur=arthur->get_point_ud_0(0);
      }
    //
    if(dry > 0)
      {
	for (int dr=0;dr < dry;dr++)
	  arthur=arthur->get_point_ud_0(3);
      }
    else if(dry < 0)
      {
	for (int dr=0;dr < -dry;dr++)
	  arthur=arthur->get_point_ud_0(2);
      }
    //
    if(drz > 0)
      {
	for (int dr=0;dr < drz;dr++)
	  arthur=arthur->get_point_ud_0(5);
      }
    else if(drz < 0)
      {
	for (int dr=0;dr < -drz;dr++)
	  arthur=arthur->get_point_ud_0(4);
      }
    return arthur;
  }
  //
  void Point::get_point_ud(vector <Point*>& point_6) const
  {
//     point_6=point_ud;
    point_6.assign(point_ud.begin(),point_ud.begin()+6);
    // if(point_6.size() < 6)
    //   point_6.resize(6);
    // std::copy(point_ud,point_ud+6,point_6.begin());
  }
  void Point::set_point_ud(vector <Point*>& point_6)
  {
//     point_ud=point_6;
    for(int ni : {0,1,2,3,4,5})
      point_ud[ni]=point_6[ni];
    // std::copy(point_6.begin(),point_6.begin()+6,point_ud);
  }
  Point* Point::get_point_up_x() const
  {
    return point_ud[1];
  }
  Point* Point::get_point_up_x_0() const
  {
    if(Point::pointer_test(point_ud[1]))
      return point_ud[1];
    *p_FILE << " up x error " << endl;
    dump();
    assert(0);
    return 0;
  }
  void Point::set_point_up_x(Point* p_point)
  {
    point_ud[1]=p_point;
  }
  Point* Point::get_point_up_y() const
  {
    return point_ud[3];
  }
  Point* Point::get_point_up_y_0() const
  {
    if(Point::pointer_test(point_ud[3]))
      return point_ud[3];
    *p_FILE << " up y error " << endl;
    dump();
    assert(0);
    return 0;
  }
  void Point::set_point_up_y(Point* p_point)
  {
    point_ud[3]=p_point;
  }
  Point* Point::get_point_up_z() const
  {
    return point_ud[5];
  }
  Point* Point::get_point_up_z_0() const
  {
    if(Point::pointer_test(point_ud[5]))
      return point_ud[5];
    *p_FILE << " up z error " << endl;
    dump();
    assert(0);
    return 0;
  }
  void Point::set_point_up_z(Point* p_point)
  {
    point_ud[5]=p_point;
  }
  Point* Point::get_point_down_x() const
  {
    return point_ud[0];
  }
  Point* Point::get_point_down_x_0() const
  {
    if(Point::pointer_test(point_ud[0]))
      return point_ud[0];
    *p_FILE << " down x error " << endl;
    dump();
    assert(0);
    return 0;
  }
  void Point::set_point_down_x(Point* p_point)
  {
    point_ud[0]=p_point;
  }
  void Point::down_from_up(Point* p_up_x,Point* p_up_y,Point* p_up_z)
  {
    if(p_up_x != 0)
      p_up_x->point_ud[0]=this;
    if(p_up_y != 0)
      p_up_y->point_ud[2]=this;
    if(p_up_z != 0)
      p_up_z->point_ud[4]=this;
  }
  void Point::down_from_up()
  {
    if(point_ud[1] != 0)
      point_ud[1]->point_ud[0]=this;
    if(point_ud[3] != 0)
      point_ud[3]->point_ud[2]=this;
    if(point_ud[5] != 0)
      point_ud[5]->point_ud[4]=this;
  }
  void Point::up_from_down()
  {
    if(point_ud[0] != 0)
      point_ud[0]->point_ud[1]=this;
    if(point_ud[2] != 0)
      point_ud[2]->point_ud[3]=this;
    if(point_ud[4] != 0)
      point_ud[4]->point_ud[5]=this;
  }
  Point* Point::get_point_down_y() const
  {
    return point_ud[2];
  }
  Point* Point::get_point_down_y_0() const
  {
    if(Point::pointer_test(point_ud[2]))
      return point_ud[2];
    *p_FILE << " down y error " << endl;
    dump();
    assert(0);
    return 0;
  }
  void Point::set_point_down_y(Point* p_point)
  {
    point_ud[2]=p_point;
  }
  Point* Point::get_point_down_z() const
  {
    return point_ud[4];
  }
  Point* Point::get_point_down_z_0() const
  {
    if(Point::pointer_test(point_ud[4]))
      return point_ud[4];
    *p_FILE << " down z error " << endl;
    dump();
    assert(0);
    return 0;
  }
  void Point::set_point_down_z(Point* p_point)
  {
    point_ud[4]=p_point;
  }
  Point* Point::get_point_ud(const int& i) const
  {
    return point_ud[i];
  }
  Point* Point::get_point_ud_0(const int& i,const int& tag) const
  {
    if(Point::pointer_test(point_ud[i]))
      return point_ud[i];
    *p_FILE << " ud error " << i << " " << tag << endl;
    dump();
    assert(0);
    return 0;
  }
  Point* Point::get_point_ud_0(const int& i) const
  {
    if(Point::pointer_test(point_ud[i]))
      return point_ud[i];
    *p_FILE << " ud error " << i << endl;
    dump();
    assert(0);
    return 0;
  }
  void Point::set_point_ud(Point* p_point,const int& i)
  {
    point_ud[i]=p_point;
  }
  Point* Point::get_point_pointer() const
  {
    if(Point::pointer_test(point_pointer))
      return point_pointer;
    *p_FILE << " point pointer error " << endl;
    dump();
    assert(point_pointer);
    return 0;
  }
  Point* Point::get_point_pointer_t() const
  {
    return point_pointer;
  }
  void Point::set_point_pointer(Point* p_point)
  {
    point_pointer=p_point;
  }
  Point* Point::get_p_daughter_point() const
  {
    bool spam=Point::pointer_test(p_daughter_point);
    if(!spam) 
      {
	*p_FILE << " daughter error " << endl;
	dump();
	assert(0);
      }
    return p_daughter_point;
  }
  void Point::set_p_daughter_point(Point* p_d_point)
  {
    p_daughter_point=p_d_point;
  }
  void Point::point_pointers_all(Point& high_point)
  {
    this->point_pointer=&high_point;

    Point* hh=high_point.get_point_ud_0(1);
    Point* p=(get_point_ud_0(1))->get_point_ud_0(1);
    p->point_pointer=hh;
    hh->p_daughter_point=p;

    hh=hh->get_point_ud_0(3);
    p=(p->get_point_ud_0(3))->get_point_ud_0(3);
    p->point_pointer=hh;
    hh->p_daughter_point=p;

    hh=hh->get_point_ud_0(0);
    p=(p->get_point_ud_0(0))->get_point_ud_0(0);
    p->point_pointer=hh;
    hh->p_daughter_point=p;

    hh=hh->get_point_ud_0(5);
    p=(p->get_point_ud_0(5))->get_point_ud_0(5);
    p->point_pointer=hh;
    hh->p_daughter_point=p;

    hh=hh->get_point_ud_0(1);
    p=(p->get_point_ud_0(1))->get_point_ud_0(1);
    p->point_pointer=hh;
    hh->p_daughter_point=p;

    hh=hh->get_point_ud_0(2);
    p=(p->get_point_ud_0(2))->get_point_ud_0(2);
    p->point_pointer=hh;
    hh->p_daughter_point=p;

    hh=hh->get_point_ud_0(0);
    p=(p->get_point_ud_0(0))->get_point_ud_0(0);
    p->point_pointer=hh;
    hh->p_daughter_point=p;
  }
  double Point::laplacian() const
  {
    return get_point_ud_0(0)->potential_point+
      get_point_ud_0(1)->potential_point+
      get_point_ud_0(2)->potential_point+
      get_point_ud_0(3)->potential_point+
      get_point_ud_0(4)->potential_point+
      get_point_ud_0(5)->potential_point-
      6.0*potential_point;
  }
//   void Point::force_shear_point_make()
//   {
//     if(force_shear_point == 0)
//       force_shear_point= new double[6];
//   }
  void Point::force_shear_point_zero()
  {
    {
      vector<double>Patsy;
      force_shear_point.swap(Patsy);
    }
    // force_shear_point.clear();
//     if(force_shear_point != 0)
//       delete [] force_shear_point;
//     force_shear_point=0;
  }
  void Point::copy_force_shear_point(Point& p0)
  {
    force_shear_point.resize(6);
    force_shear_point=p0.force_shear_point;
//     std::copy(p0.force_shear_point,p0.force_shear_point+6,force_shear_point);
  }
  void Point::copy_density_point(Point& p0,Point& p1)
  {
    density_point=0.5*(p0.density_point+p1.density_point);
  }
  void Point::copy_density_point(Point& p0,Point& p1,Point& p2,Point& p3)
  {
    density_point=0.25*(p0.density_point+p1.density_point+p2.density_point+p3.density_point);
  }
  void Point::copy_potential_point(Point& p0)
  {
    potential_point=p0.potential_point;
  }
  void Point::copy_force_shear_point_6()
  {
    force_shear_point.resize(6);
    for(unsigned int ni=0;ni<6;ni++)
      force_shear_point[ni]=(get_point_ud_0(0)->force_shear_point[ni]+get_point_ud_0(1)->force_shear_point[ni]+get_point_ud_0(2)->force_shear_point[ni]+
			     get_point_ud_0(3)->force_shear_point[ni]+get_point_ud_0(4)->force_shear_point[ni]+get_point_ud_0(5)->force_shear_point[ni])/6.0;
  }
  void Point::copy_force_shear_point_4(vector <int>& witch)
  {
    force_shear_point.resize(6);
    for(unsigned int ni=0;ni<6;ni++)
      force_shear_point[ni]=(get_point_ud_0(witch[0])->force_shear_point[ni]+get_point_ud_0(witch[1])->force_shear_point[ni]+
			     get_point_ud_0(witch[2])->force_shear_point[ni]+get_point_ud_0(witch[3])->force_shear_point[ni])*0.25;
  }
  void Point::copy_force_shear_point_2(vector <int>& witch)
  {
    force_shear_point.resize(6);
    for(unsigned int ni=0;ni<6;ni++)
      force_shear_point[ni]=(get_point_ud_0(witch[0])->force_shear_point[ni]+get_point_ud_0(witch[1])->force_shear_point[ni])*0.5;
  }
  void Point::copy_force_shear_point_1()
  {
    force_shear_point.resize(6);
    if(point_pointer != 0)
      {
	force_shear_point=point_pointer->force_shear_point;
// 	std::copy(point_pointer->force_shear_point,point_pointer->force_shear_point+6,force_shear_point);
	return;
      }
    dump();
    assert(point_pointer);
  }
  void Point::copy_force_point_6()
  {
    for(int ni=0;ni<3;ni++)
      force_point[ni]=(get_point_ud_0(0)->force_point[ni]+get_point_ud_0(1)->force_point[ni]+get_point_ud_0(2)->force_point[ni]+
		       get_point_ud_0(3)->force_point[ni]+get_point_ud_0(4)->force_point[ni]+get_point_ud_0(5)->force_point[ni])/6.0;
  }
  void Point::copy_force_point_4(vector <int>& witch)
  {
    for(int ni=0;ni<3;ni++)
      force_point[ni]=(get_point_ud_0(witch[0])->force_point[ni]+get_point_ud_0(witch[1])->force_point[ni]+
		       get_point_ud_0(witch[2])->force_point[ni]+get_point_ud_0(witch[3])->force_point[ni])*0.25;
  }
  void Point::copy_force_point_2(vector <int>& witch)
  {
    for(int ni=0;ni<3;ni++)
      force_point[ni]=(get_point_ud_0(witch[0])->force_point[ni]+get_point_ud_0(witch[1])->force_point[ni])*0.5;
  }
  void Point::copy_force_point_1()
  {
    if(point_pointer != 0)
      {
	for(int ni=0;ni<3;ni++)
	  force_point[ni]=point_pointer->force_point[ni];
	return;
      }
    dump();
    assert(point_pointer);
  }
  void Point::copy_potential_point_6()
  {
    potential_point=(get_point_ud_0(0)->potential_point+get_point_ud_0(1)->potential_point+get_point_ud_0(2)->potential_point+
		     get_point_ud_0(3)->potential_point+get_point_ud_0(4)->potential_point+get_point_ud_0(5)->potential_point)/6.0;
  }
  void Point::copy_potential_point_4(vector <int>& witch)
  {
    potential_point=(get_point_ud_0(witch[0])->potential_point+get_point_ud_0(witch[1])->potential_point+
		     get_point_ud_0(witch[2])->potential_point+get_point_ud_0(witch[3])->potential_point)*0.25;
  }
  void Point::copy_potential_point_2(vector <int>& witch)
  {
    potential_point=(get_point_ud_0(witch[0])->potential_point+get_point_ud_0(witch[1])->potential_point)*0.5;
  }
  void Point::copy_potential_point_1()
  {
    if(point_pointer != 0)
      {
	potential_point=point_pointer->potential_point;
	return;
      }
    dump();
    assert(point_pointer);
  }
  void Point::diff_pot(const double& conv)
  {
    force_point[0]=(get_point_ud_0(0)->potential_point-get_point_ud_0(1)->potential_point)*conv;
    force_point[1]=(get_point_ud_0(2)->potential_point-get_point_ud_0(3)->potential_point)*conv;
    force_point[2]=(get_point_ud_0(4)->potential_point-get_point_ud_0(5)->potential_point)*conv;
  }
  void Point::diff_pot_careful(const double& conv)
  {
    if(point_ud[0] != 0 && point_ud[1] != 0 && point_ud[2] != 0 && point_ud[3] != 0 && point_ud[4] != 0 && point_ud[5] != 0)
      diff_pot(conv);
  }
  void Point::diff_force(const double& conv)
  {
    // (fx,x)(fx,y)(fx,z)(fy,y)(fy,z)(fz,z)
    force_shear_point.resize(6);
    force_shear_point[0]=-(get_point_ud_0(0)->force_point[0]-get_point_ud_0(1)->force_point[0])*conv;
    force_shear_point[1]=-(get_point_ud_0(2)->force_point[0]-get_point_ud_0(3)->force_point[0])*conv;
    force_shear_point[2]=-(get_point_ud_0(4)->force_point[0]-get_point_ud_0(5)->force_point[0])*conv;
    force_shear_point[3]=-(get_point_ud_0(2)->force_point[1]-get_point_ud_0(3)->force_point[1])*conv;
    force_shear_point[4]=-(get_point_ud_0(4)->force_point[1]-get_point_ud_0(5)->force_point[1])*conv;
    force_shear_point[5]=-(get_point_ud_0(4)->force_point[2]-get_point_ud_0(5)->force_point[2])*conv;
  }
  void Point::div_force(const double& conv)
  {
    // (fx,x)(fx,y)(fx,z)(fy,y)(fy,z)(fz,z)
    force_shear_point.resize(1);
//     if(force_shear_point != 0)
//       delete [] force_shear_point;
//     force_shear_point= new double[1];
    force_shear_point[0]=
      (-(get_point_ud_0(0)->force_point[0]-get_point_ud_0(1)->force_point[0])
       -(get_point_ud_0(2)->force_point[0]-get_point_ud_0(3)->force_point[0])
       -(get_point_ud_0(4)->force_point[0]-get_point_ud_0(5)->force_point[0]))*conv;
  }
  double Point::get_div_force() const
  {
    return force_shear_point[0];
  }
  void Point::get_shear(vector <double>& s) const
  {
    s=force_shear_point;
//     if(s.size() < 6)
//       s.resize(6);
//     std::copy(force_shear_point,force_shear_point+6,s.begin());
  }
  void Point::add_dens_at_point(const double& d)
  {
    density_point+=d;
  }
  template <class T> void Point::add_density_at_points(vector <T>& d)
  {
    density_point+=d[0];
    Point* p=get_point_ud_0(1);
    p->density_point+=d[1];
    p=p->get_point_ud_0(3);
    p->density_point+=d[3];
    p=p->get_point_ud_0(0);
    p->density_point+=d[2];
    p=p->get_point_ud_0(5);
    p->density_point+=d[6];
    p=p->get_point_ud_0(1);
    p->density_point+=d[7];
    p=p->get_point_ud_0(2);
    p->density_point+=d[5];
    p=p->get_point_ud_0(0);
    p->density_point+=d[4];
  }
  template void Point::add_density_at_points(vector <double>& d);
  void Point::subtract_dens_at_point(const double& d)
  {
    density_point-=d;
    //      *p_FILE << " subtract " << pos_point[0] << " " << pos_point[1] << " " << pos_point[2] << " " << inside << "\n";
  }
  void Point::scale_density_point(const double& s)
  {
    density_point*=s;
  }
  void Point::dens_from_mother()
  {
    density_point=point_pointer->density_point;
  }
  void Point::set_density_point(const double& dens)
  {
    density_point=dens;
  }
  void Point::get_field_values(vector <double>& pott,vector <double>& fx,vector <double>& fy,vector <double>& fz) const
  {
    pott[0]=potential_point;
    fx[0]=force_point[0];
    fy[0]=force_point[1];
    fz[0]=force_point[2];
    Point* p=get_point_ud_0(1);
    pott[1]=p->potential_point;
    fx[1]=p->force_point[0];
    fy[1]=p->force_point[1];
    fz[1]=p->force_point[2];
    p=p->get_point_ud_0(3);
    pott[3]=p->potential_point;
    fx[3]=p->force_point[0];
    fy[3]=p->force_point[1];
    fz[3]=p->force_point[2];
    p=p->get_point_ud_0(0);
    pott[2]=p->potential_point;
    fx[2]=p->force_point[0];
    fy[2]=p->force_point[1];
    fz[2]=p->force_point[2];
    p=p->get_point_ud_0(5);
    pott[6]=p->potential_point;
    fx[6]=p->force_point[0];
    fy[6]=p->force_point[1];
    fz[6]=p->force_point[2];
    p=p->get_point_ud_0(1);
    pott[7]=p->potential_point;
    fx[7]=p->force_point[0];
    fy[7]=p->force_point[1];
    fz[7]=p->force_point[2];
    p=p->get_point_ud_0(2);
    pott[5]=p->potential_point;
    fx[5]=p->force_point[0];
    fy[5]=p->force_point[1];
    fz[5]=p->force_point[2];
    p=p->get_point_ud_0(0);
    pott[4]=p->potential_point;
    fx[4]=p->force_point[0];
    fy[4]=p->force_point[1];
    fz[4]=p->force_point[2];
  }
  void Point::get_field_shear_values(vector <double>& fxx,vector <double>& fxy,vector <double>& fxz,
			      vector <double>& fyy,vector <double>& fyz,vector <double>& fzz) const
  {
    fxx[0]=force_shear_point[0];
    fxy[0]=force_shear_point[1];
    fxz[0]=force_shear_point[2];
    fyy[0]=force_shear_point[3];
    fyz[0]=force_shear_point[4];
    fzz[0]=force_shear_point[5];
    Point* p=get_point_ud_0(1);
    fxx[1]=p->force_shear_point[0];
    fxy[1]=p->force_shear_point[1];
    fxz[1]=p->force_shear_point[2];
    fyy[1]=p->force_shear_point[3];
    fyz[1]=p->force_shear_point[4];
    fzz[1]=p->force_shear_point[5];
    p=p->get_point_ud_0(3);
    fxx[3]=p->force_shear_point[0];
    fxy[3]=p->force_shear_point[1];
    fxz[3]=p->force_shear_point[2];
    fyy[3]=p->force_shear_point[3];
    fyz[3]=p->force_shear_point[4];
    fzz[3]=p->force_shear_point[5];
    p=p->get_point_ud_0(0);
    fxx[2]=p->force_shear_point[0];
    fxy[2]=p->force_shear_point[1];
    fxz[2]=p->force_shear_point[2];
    fyy[2]=p->force_shear_point[3];
    fyz[2]=p->force_shear_point[4];
    fzz[2]=p->force_shear_point[5];
    p=p->get_point_ud_0(5);
    fxx[6]=p->force_shear_point[0];
    fxy[6]=p->force_shear_point[1];
    fxz[6]=p->force_shear_point[2];
    fyy[6]=p->force_shear_point[3];
    fyz[6]=p->force_shear_point[4];
    fzz[6]=p->force_shear_point[5];
    p=p->get_point_ud_0(1);
    fxx[7]=p->force_shear_point[0];
    fxy[7]=p->force_shear_point[1];
    fxz[7]=p->force_shear_point[2];
    fyy[7]=p->force_shear_point[3];
    fyz[7]=p->force_shear_point[4];
    fzz[7]=p->force_shear_point[5];
    p=p->get_point_ud_0(2);
    fxx[5]=p->force_shear_point[0];
    fxy[5]=p->force_shear_point[1];
    fxz[5]=p->force_shear_point[2];
    fyy[5]=p->force_shear_point[3];
    fyz[5]=p->force_shear_point[4];
    fzz[5]=p->force_shear_point[5];
    p=p->get_point_ud_0(0);
    fxx[4]=p->force_shear_point[0];
    fxy[4]=p->force_shear_point[1];
    fxz[4]=p->force_shear_point[2];
    fyy[4]=p->force_shear_point[3];
    fyz[4]=p->force_shear_point[4];
    fzz[4]=p->force_shear_point[5];
  }
  void Point::clean_shear()
  {
    vector<double>King_Arthur;
    force_shear_point.swap(King_Arthur);
  }
  void Point::get_field_values(vector <double>& pott) const
  {
    pott[0]=potential_point;
    Point* p=get_point_ud_0(1);
    pott[1]=p->potential_point;
    p=p->get_point_ud_0(3);
    pott[3]=p->potential_point;
    p=p->get_point_ud_0(0);
    pott[2]=p->potential_point;
    p=p->get_point_ud_0(5);
    pott[6]=p->potential_point;
    p=p->get_point_ud_0(1);
    pott[7]=p->potential_point;
    p=p->get_point_ud_0(2);
    pott[5]=p->potential_point;
    p=p->get_point_ud_0(0);
    pott[4]=p->potential_point;
  }
  void Point::get_density_points(vector <double>& dens)
  {
    dens[0]=density_point;
    Point* p=get_point_ud_0(1);
    dens[1]=p->density_point;
    p=p->get_point_ud_0(3);
    dens[3]=p->density_point;
    p=p->get_point_ud_0(0);
    dens[2]=p->density_point;
    p=p->get_point_ud_0(5);
    dens[6]=p->density_point;
    p=p->get_point_ud_0(1);
    dens[7]=p->density_point;
    p=p->get_point_ud_0(2);
    dens[5]=p->density_point;
    p=p->get_point_ud_0(0);
    dens[4]=p->density_point;
  }
  void Point::dumpy() const
  {
    *p_FILE << "dumpy " << this << " ";
    *p_FILE << real_pointer << " ";
    *p_FILE << pos_point[0] << " ";
    *p_FILE << pos_point[1] << " ";
    *p_FILE << pos_point[2] << "\n";
  }
  void Point::dumpd() const
  {
    *p_FILE << "dumpd " << this << " ";
    *p_FILE << real_pointer << " ";
    *p_FILE << pos_point[0] << " ";
    *p_FILE << pos_point[1] << " ";
    *p_FILE << pos_point[2] << " ";
    *p_FILE << inside << " ";
    *p_FILE << edge_point << " ";
    *p_FILE << buffer_point << " ";
    *p_FILE << passive_point << " ";
    *p_FILE << density_point << "\n";
  }
  void Point::dumpp() const
  {
    *p_FILE << "dumpp " << this << " ";
    *p_FILE << real_pointer << " ";
    *p_FILE << inside << " ";
    *p_FILE << edge_point << " ";
    *p_FILE << buffer_point << " ";
    *p_FILE << passive_point << " ";
    *p_FILE << pos_point[0] << " ";
    *p_FILE << pos_point[1] << " ";
    *p_FILE << pos_point[2] << " ";
    *p_FILE << setprecision(6) << density_point << " ";
    *p_FILE << potential_point << "\n";
  }
  void Point::dumpp(ofstream& FF) const
  {
    FF << "dumpp " << this << " ";
    FF << real_pointer << " ";
    FF << inside << " ";
    FF << edge_point << " ";
    FF << buffer_point << " ";
    FF << passive_point << " ";
    FF << pos_point[0] << " ";
    FF << pos_point[1] << " ";
    FF << pos_point[2] << " ";
    FF << setprecision(6) << density_point << " ";
    FF << potential_point << "\n";
  }
  void Point::dumppf() const
  {
    *p_FILE << "dumppf " << this << " ";
    *p_FILE << inside << " ";
    *p_FILE << edge_point << " ";
    *p_FILE << buffer_point << " ";
    *p_FILE << passive_point << " ";
    *p_FILE << real_pointer << " ";
    *p_FILE << pos_point[0] << " ";
    *p_FILE << pos_point[1] << " ";
    *p_FILE << pos_point[2] << " ";
    *p_FILE << setprecision(6) << density_point << " ";
    *p_FILE << potential_point << " ";
    *p_FILE << force_point[0] << " ";
    *p_FILE << force_point[1] << " ";
    *p_FILE << force_point[2] << "\n";
  }
  void Point::dump() const
  {
    *p_FILE << "\t" << this << "\t";
    *p_FILE << pos_point[0] << "\t";
    *p_FILE << pos_point[1] << "\t";
    *p_FILE << pos_point[2] << "\t";
    *p_FILE << "dump point " << " ";
    *p_FILE << inside << " ";
    *p_FILE << edge_point << " ";
    *p_FILE << buffer_point << " ";
    *p_FILE << passive_point << " ";
    *p_FILE << real_pointer << "\n";
    //
    for(int ni=0;ni<6;ni++)
      {
	*p_FILE << ni << "\t" << point_ud[ni] << "\t";
	if(point_ud[ni] != 0) 
	  *p_FILE << point_ud[ni]->pos_point[0] << "\t" << point_ud[ni]->pos_point[1] << "\t" << point_ud[ni]->pos_point[2];
	*p_FILE << "\n";
      }
  } 
  void Point::get_field_values(double& pot,double& fx,double& fy,double& fz) const
  {
    pot=potential_point;
    fx=force_point[0];
    fy=force_point[1];
    fz=force_point[2];
  }
  void Point::scale_pot_forces(const double& scaling)
  {
    potential_point*=scaling;
    force_point[0]*=scaling;
    force_point[1]*=scaling;
    force_point[2]*=scaling;
  }
  void Point::get_deltas(vector <double>& pos,double& d_x, double& d_y, double& d_z,const double& scale,const double& d_inv) const
  {
    d_x=(pos[0]*scale-pos_point[0])*d_inv;
    d_y=(pos[1]*scale-pos_point[1])*d_inv;
    d_z=(pos[2]*scale-pos_point[2])*d_inv;
    if(abs(d_x-0.5) > 0.5 || abs(d_y-0.5) > 0.5 || abs(d_z-0.5) > 0.5)
      {
	*p_FILE << " this " << this << "\n";
	*p_FILE << "deltax " << pos[0] << " " << pos_point[0] << " " << scale << " " << d_inv << "\n";
	*p_FILE << "deltay " << pos[1] << " " << pos_point[1] << " " << scale << " " << d_inv << "\n";
	*p_FILE << "deltaz " << pos[2] << " " << pos_point[2] << " " << scale << " " << d_inv << "\n";
	*p_FILE << "dxyz " << d_x << " " << d_y << " " << d_z << endl;
      }
//     assert(abs(d_x-0.5) <= 0.5);
//     assert(abs(d_y-0.5) <= 0.5);
//     assert(abs(d_z-0.5) <= 0.5);
  }
  void Point::get_potss_denss(const int* dprl,const double& g_c,double* potss,double* denss)
  {
    int n=0;
    Point* p_point_z=this;
    for(int nz=1;nz <= dprl[2];++nz)
      {
	Point* p_point_y=p_point_z;
	for(int ny=1;ny <= dprl[1];++ny)
	  {
	    Point* p_point_x=p_point_y; 
	    for(int nx=1;nx <= dprl[0];++nx)
	      {
		potss[n]=p_point_x->get_potential_point();
		denss[n]=g_c*p_point_x->get_density_point();
		if(nx < dprl[0])
		  p_point_x=p_point_x->get_point_up_x_0();
		n++;
	      }
	    if(ny < dprl[1])
	      p_point_y=p_point_y->get_point_up_y_0();
	  }
	if(nz < dprl[2])
	  p_point_z=p_point_z->get_point_up_z_0();
      }
  }
  void Point::set_potss(const int* dprl,const double* potss)
  {
    int n=0;
    Point* p_point_z=this;
    for(int nz=1;nz <= dprl[2];++nz)
      {
	Point* p_point_y=p_point_z;
	for(int ny=1;ny <= dprl[1];++ny)
	  {
	    Point* p_point_x=p_point_y; 
	    for(int nx=1;nx <= dprl[0];++nx)
	      {
		p_point_x->set_potential_point(potss[n]);
		if(nx < dprl[0])
		  p_point_x=p_point_x->get_point_up_x_0();
		n++;
	      }
	    if(ny < dprl[1])
	      p_point_y=p_point_y->get_point_up_y_0();
	  }
	if(nz < dprl[2])
	  p_point_z=p_point_z->get_point_up_z_0();
      }
  }
  Point* Point::neighbor_high(const int& what) const
  {
    if(point_ud[what] != 0 && point_ud[what]->get_it_is_high()) return point_ud[what];
    return 0;
  }
  Point* Point::neighbor_high(const int& what1,const int& what2) const
  {
    Point* ref=neighbor_high(what1);
    if(ref != 0) return ref->neighbor_high(what2);
    ref=neighbor_high(what2);
    if(ref != 0) return ref->neighbor_high(what1);
    return 0;
  }
  Point* Point::neighbor_high(const int& what1,const int& what2,const int& what3) const
  {
    Point* ref=neighbor_high(what1,what2);
    if(ref != 0) return ref->neighbor_high(what3);
    ref=neighbor_high(what2,what3);
    if(ref != 0) return ref->neighbor_high(what1);
    ref=neighbor_high(what1,what3);
    if(ref != 0) return ref->neighbor_high(what2);
    return 0;
  }
  Point* Point::neighbor(const int& what) const
  {
    return point_ud[what];
  }
  Point* Point::neighbor(const int& what1,const int& what2) const
  {
    Point* ref=neighbor(what1);
    if(ref != 0) return ref->neighbor(what2);
    ref=neighbor(what2);
    if(ref != 0) return ref->neighbor(what1);
    return 0;
  }
  Point* Point::neighbor(const int& what1,const int& what2,const int& what3) const
  {
    Point* ref=neighbor(what1,what2);
    if(ref != 0) return ref->neighbor(what3);
    ref=neighbor(what2,what3);
    if(ref != 0) return ref->neighbor(what1);
    ref=neighbor(what1,what3);
    if(ref != 0) return ref->neighbor(what2);
    return 0;
  }
  double Point::sum_inside_weights(vector <double>& w)
  {
    double sum=1.0e-30;
    Point* p=this;
    for(int i=0;i<8;i++)
      {
	if(p->inside)
	  sum+=w[Point::sequence[i]];
	p=p->point_ud[Point::updown[i]];
      }
    return sum;
  }
}
