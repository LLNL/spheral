#ifndef _Point_Defined_
#define _Point_Defined_
namespace FractalSpace
{
  class Point{
    int real_pointer;
    int point_to_number;
    int high_number;
    Point* point_pointer;
    Point* p_daughter_point;
    Group* p_in_group;
    Group* p_in_high_group;
    bool inside;
    bool it_is_high;
    //
    bool buffer_point;
    bool surface_point;
    int rank_point;
    //
    vector <bool> eureka_adj;
    vector <bool> eureka_dau;
    double potential_point;
    double density_point;
    vector <int> pos_point;
    vector <Point*> point_ud;
    vector <double>force_point;
    vector <double>force_shear_point;
  public:
    vector <Particle*> list_particles;
    vector <Particle*> list_other_particles;
    static Point* nothing;
    static int number_points;
    static bool calc_candidates;
    static int order[8][7];
    static vector <int> phl;
    static vector <int> dupes;
    static vector <int> corner_a;
    static vector <int> corner_b;
    static vector <vector <int> >nextt;
    static vector <bool> left;
    static vector <bool> corner;
    static vector <bool> edge;
    static vector <bool> face;
    static vector <bool> center;
    static vector <int> sequence;
    static vector <int> updown;
    static vector < vector <int> > cefc;
    Point():
      real_pointer(-1),
      point_to_number(-1),
      high_number(-1),
      point_pointer(0),
      p_daughter_point(0),
      p_in_group(0),
      p_in_high_group(0),
      inside(false),
      it_is_high(false),
      buffer_point(false),
      surface_point(false),
      rank_point(0),
      potential_point(0.0),
      density_point(0.0),
      pos_point(3,-1),
      point_ud(6),
      force_point(3,0.0)
    {
      for(int ii=0;ii<6;ii++)
	point_ud[ii]=0;
      number_points++;
    };
    Point(const int& r_p,const bool& b_p,const bool& s_p):
      real_pointer(-1),
      point_to_number(-1),
      high_number(-1),
      point_pointer(0),
      p_daughter_point(0),
      p_in_group(0),
      p_in_high_group(0),
      inside(false),
      it_is_high(false),
      potential_point(0.0),
      density_point(0.0),
      pos_point(3,-1),
      point_ud(6),
      force_point(3,0.0)
    {
      rank_point=r_p;
      buffer_point=b_p;
      surface_point=s_p;
      for(int ii=0;ii<6;ii++)
	point_ud[ii]=0;
      number_points++;
    };
    ~Point()
    {    
      number_points--;
    };
    /*
    void clean_point()
    {
      real_pointer=-1;
      point_to_number=-1;
      high_number=-1;
      point_pointer=0;
      p_daughter_point=0;
      p_in_group=0;
      p_in_high_group=0;
      inside=false;
      it_is_high=false;
      eureka_adj.assign(27,false);
      eureka_dau.assign(27,false);
      potential_point=0.0;
      density_point=0.0;
      pos_point.assign(3,-1);
      force_point.assign(3,0.0);
      force_shear_point.clear();
      point_ud.assign(6,nothing);
      list_particles.clear();
      list_other_particles.clear();
    };
    */
    void set_rank_point(const int& r_p)
    {
      rank_point=r_p;
    }
    int get_rank_point()
    {
      return rank_point;
    }
    void set_buffer_surface_point(const bool& b_p,const bool& s_p)
    {
      buffer_point=b_p;
      surface_point=s_p;
    }
    void get_buffer_surface_point(bool& b_p,bool& s_p)
    {
      b_p=buffer_point;
      s_p=surface_point;
    }
    void set_buffer_point(const bool& b_p)
    {
      buffer_point=b_p;
    }
    bool get_buffer_point()
    {
      return buffer_point;
    }
    void set_surface_point(const bool& s_p)
    {
      surface_point=s_p;
    }
    bool get_surface_point()
    {
      return surface_point;
    }
    void set_eureka_adj(const vector <bool>& eu)
    {
      eureka_adj=eu;
    }
    void get_eureka_adj(vector <bool>& eu)
    {
      eu=eureka_adj;
    }
    void set_eureka_dau(const vector <bool>& eu)
    {
      eureka_dau=eu;
    }
    void get_eureka_dau(vector <bool>& eu)
    {
      eu=eureka_dau;
    }
    bool get_eureka_dau(const int& ni)
    {
      return eureka_dau[ni];
    }
    void set_p_in_group(Group* p_g) 
    {
      p_in_group=p_g;
    }
    Group* get_p_in_group()
    {
      return p_in_group;
    }
    void set_p_in_high_group(Group* p_g) 
    {
      p_in_high_group=p_g;
    }
    Group* get_p_in_high_group()
    {
      return p_in_high_group;
    }
    void get_pos_point(vector <int>& pos)
    {
      pos=pos_point;
    }
    void set_pos_point(const vector <int>& pos)
    {
      pos_point=pos;
    }
    void set_pos_point(const int& x,const int& y,const int& z)
    {
      pos_point[0]=x;
      pos_point[1]=y;
      pos_point[2]=z;
      //      cout << "point " << this << " " << x << " "  << y << " "  << z << endl;
    }
    void get_pos_point(int& x,int& y,int& z)
    {
      x=pos_point[0];
      y=pos_point[1];
      z=pos_point[2];
    }
    int get_pos_point_x()
    {
      return pos_point[0];
    }
    void set_pos_point_x(const int& i)
    {
      pos_point[0]=i;
    }
    int get_pos_point_y()
    {
      return pos_point[1];
    }
    void set_pos_point_y(const int& i)
    {
      pos_point[1]=i;
    }
    int get_pos_point_z()
    {
      return pos_point[2];
    }
    void set_pos_point_z(const int& i)
    {
      pos_point[2]=i;
    }
    int get_pos_point(const int& i)
    {
      return pos_point[i];
    }
    int get_real_pointer()
    {
      return real_pointer;
    }
    void set_real_pointer(const int& i)
    {
      real_pointer=i;
    }
    void set_point_to_number(const int& i)
    {
      point_to_number=i;
    }
    int get_point_to_number()
    {
      return point_to_number;
    }
    void set_high_number(const int& i)
    {
      high_number=i;
    }
    int get_high_number()
    {
      return high_number;
    }
    bool get_inside()
    {
      return  inside;
    }
    void set_inside(const bool& ins)
    {
      inside=ins;
    }
    bool get_it_is_high()
    {
      return it_is_high;
    }
    void set_it_is_high(const bool& value)
    {
      it_is_high=value;
    }
    void set_inside_high()
    {
      it_is_high=it_is_high && inside; 
    }
    double get_potential_point()
    {
      return potential_point;
    }
    void set_potential_point(const double& d)
    {
      potential_point=d;
    }
    void add_potential_point(const double& d)
    {
      potential_point+=d;
    }
    void add_force_point_x(const double& d)
    {
      force_point[0]+=d;
    }
    void add_force_point_y(const double& d)
    {
      force_point[1]+=d;
    }
    void add_force_point_z(const double& d)
    {
      force_point[2]+=d;
    }
    double get_density_point()
    {
      return density_point;
    }
    void get_force_point(vector <double>& force)
    {
      force=force_point;
    }
    void set_force_point(vector <double>& force)
    {
      force_point=force;
    }
    double get_force_point_x()
    {
      return force_point[0];
    }
    //
    double get_force_point_y()
    {
      return force_point[1];
    }
    //
    double get_force_point_z()
    {
      return force_point[2];
    }
    void set_force_point_x(const double& f)
    {
      force_point[0]=f;
    }
    void set_force_point_y(const double& f)
    {
      force_point[1]=f;
    }
    void set_force_point_z(const double& f)
    {
      force_point[2]=f;
    }
    Point* move_adj(const int& ra,const int& rb)
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
    Point* move_rp(const int& r)
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
    void get_point_ud(vector <Point*>& point_6)
    {
      point_6=point_ud;
    }
    void set_point_ud(vector <Point*>& point_6)
    {
      point_ud=point_6;
    }
    Point* get_point_up_x()
    {
      return point_ud[1];
    }
    Point* get_point_up_x_0()
    {
      if(Point::pointer_test(point_ud[1]))
	return point_ud[1];
      dump();
      assert(0);
      return 0;
    }
    void set_point_up_x(Point* p_point)
    {
      point_ud[1]=p_point;
    }
    Point* get_point_up_y()
    {
      return point_ud[3];
    }
    Point* get_point_up_y_0()
    {
      if(Point::pointer_test(point_ud[3]))
	return point_ud[3];
      dump();
      assert(0);
      return 0;
    }
    void set_point_up_y(Point* p_point)
    {
      point_ud[3]=p_point;
    }
    Point* get_point_up_z()
    {
      return point_ud[5];
    }
    Point* get_point_up_z_0()
    {
      if(Point::pointer_test(point_ud[5]))
	return point_ud[5];
      dump();
      assert(0);
      return 0;
    }
    void set_point_up_z(Point* p_point)
    {
      point_ud[5]=p_point;
    }
    Point* get_point_down_x()
    {
      return point_ud[0];
    }
    Point* get_point_down_x_0()
    {
      if(Point::pointer_test(point_ud[0]))
	return point_ud[0];
      dump();
      assert(0);
      return 0;
    }
    void set_point_down_x(Point* p_point)
    {
      point_ud[0]=p_point;
    }
    void down_from_up(Point* p_up_x,Point* p_up_y,Point* p_up_z)
    {
      if(p_up_x != 0)
	p_up_x->point_ud[0]=this;
      if(p_up_y != 0)
	p_up_y->point_ud[2]=this;
      if(p_up_z != 0)
	p_up_z->point_ud[4]=this;
    }
    void down_from_up()
    {
      if(point_ud[1] != 0)
	point_ud[1]->point_ud[0]=this;
      if(point_ud[3] != 0)
	point_ud[3]->point_ud[2]=this;
      if(point_ud[5] != 0)
	point_ud[5]->point_ud[4]=this;
    }
    Point* get_point_down_y()
    {
      return point_ud[2];
    }
    Point* get_point_down_y_0()
    {
      if(Point::pointer_test(point_ud[2]))
	return point_ud[2];
      dump();
      assert(0);
      return 0;
    }
    void set_point_down_y(Point* p_point)
    {
      point_ud[2]=p_point;
    }
    Point* get_point_down_z()
    {
      return point_ud[4];
    }
    Point* get_point_down_z_0()
    {
      if(Point::pointer_test(point_ud[4]))
	return point_ud[4];
      dump();
      assert(0);
      return 0;
    }
    void set_point_down_z(Point* p_point)
    {
      point_ud[4]=p_point;
    }
    Point* get_point_ud(const int& i)
    {
      return point_ud[i];
    }
    Point* get_point_ud_0(const int& i)
    {
      if(Point::pointer_test(point_ud[i]))
	return point_ud[i];
      dump();
      assert(0);
      return 0;
    }
    void set_point_ud(Point* p_point,const int& i)
    {
      point_ud[i]=p_point;
    }
    Point* get_point_pointer()
    {
      if(Point::pointer_test(point_pointer))
	return point_pointer;
      dump();
      assert(point_pointer);
      return 0;
    }
    Point* get_point_pointer_t()
    {
      return point_pointer;
    }
    void set_point_pointer(Point* p_point)
    {
      point_pointer=p_point;
    }
    Point* get_p_daughter_point()
    {
      bool spam=Point::pointer_test(p_daughter_point);
      if(!spam) 
	{
	  dump();
	  assert(0);
	}
      return p_daughter_point;
    }
    void set_p_daughter_point(Point* p_d_point)
    {
      p_daughter_point=p_d_point;
    }
    void point_pointers_all(Point& high_point)
    {
      this->point_pointer=&high_point;

      Point* hh=high_point.get_point_ud_0(1);
      Point* p=(get_point_ud_0(1))->get_point_ud_0(1);
      p->point_pointer=hh;

      hh=hh->get_point_ud_0(3);
      p=(p->get_point_ud_0(3))->get_point_ud_0(3);
      p->point_pointer=hh;

      hh=hh->get_point_ud_0(0);
      p=(p->get_point_ud_0(0))->get_point_ud_0(0);
      p->point_pointer=hh;

      hh=hh->get_point_ud_0(5);
      p=(p->get_point_ud_0(5))->get_point_ud_0(5);
      p->point_pointer=hh;

      hh=hh->get_point_ud_0(1);
      p=(p->get_point_ud_0(1))->get_point_ud_0(1);
      p->point_pointer=hh;

      hh=hh->get_point_ud_0(2);
      p=(p->get_point_ud_0(2))->get_point_ud_0(2);
      p->point_pointer=hh;

      hh=hh->get_point_ud_0(0);
      p=(p->get_point_ud_0(0))->get_point_ud_0(0);
      p->point_pointer=hh;
    }
    double laplacian()
    {
      return get_point_ud_0(0)->potential_point+
	get_point_ud_0(1)->potential_point+
	get_point_ud_0(2)->potential_point+
	get_point_ud_0(3)->potential_point+
	get_point_ud_0(4)->potential_point+
	get_point_ud_0(5)->potential_point-
	6.0*potential_point;
    }
    void force_shear_point_zero()
    {
      force_shear_point.clear();
    }
    void copy_force_shear_point(Point& p0)
    {
      force_shear_point.resize(6);
      force_shear_point=p0.force_shear_point;
    }
    void copy_density_point(Point& p0,Point& p1)
    {
      density_point=0.5*(p0.density_point+p1.density_point);
    }
    void copy_density_point(Point& p0,Point& p1,Point& p2,Point& p3)
    {
      density_point=0.25*(p0.density_point+p1.density_point+p2.density_point+p3.density_point);
    }
    void copy_potential_point(Point& p0)
    {
      potential_point=p0.potential_point;
    }
    void copy_force_shear_point_6()
    {
      force_shear_point.resize(6);
      for(unsigned int ni=0;ni<force_shear_point.size();ni++)
	force_shear_point[ni]=(point_ud[0]->force_shear_point[ni]+point_ud[1]->force_shear_point[ni]+point_ud[2]->force_shear_point[ni]+
		       point_ud[3]->force_shear_point[ni]+point_ud[4]->force_shear_point[ni]+point_ud[5]->force_shear_point[ni])/6.0;
    }
    void copy_force_shear_point_4(vector <int>& witch)
    {
      force_shear_point.resize(6);
      for(unsigned int ni=0;ni<force_shear_point.size();ni++)
	force_shear_point[ni]=(point_ud[witch[0]]->force_shear_point[ni]+point_ud[witch[1]]->force_shear_point[ni]+
		       point_ud[witch[2]]->force_shear_point[ni]+point_ud[witch[3]]->force_shear_point[ni])*0.25;
    }
    void copy_force_shear_point_2(vector <int>& witch)
    {
      force_shear_point.resize(6);
      for(unsigned int ni=0;ni<force_shear_point.size();ni++)
	force_shear_point[ni]=(point_ud[witch[0]]->force_shear_point[ni]+point_ud[witch[1]]->force_shear_point[ni])*0.5;
    }
    void copy_force_shear_point_1()
    {
      force_shear_point.resize(6);
      if(point_pointer != 0)
	{
	  for(unsigned int ni=0;ni<force_shear_point.size();ni++)
	    force_shear_point[ni]=point_pointer->force_shear_point[ni];
	  return;
	}
      dump();
      assert(point_pointer);
    }
    void copy_force_point_6()
    {
      for(int ni=0;ni<3;ni++)
	force_point[ni]=(point_ud[0]->force_point[ni]+point_ud[1]->force_point[ni]+point_ud[2]->force_point[ni]+
		       point_ud[3]->force_point[ni]+point_ud[4]->force_point[ni]+point_ud[5]->force_point[ni])/6.0;
    }
    void copy_force_point_4(vector <int>& witch)
    {
      for(int ni=0;ni<3;ni++)
	force_point[ni]=(point_ud[witch[0]]->force_point[ni]+point_ud[witch[1]]->force_point[ni]+
		       point_ud[witch[2]]->force_point[ni]+point_ud[witch[3]]->force_point[ni])*0.25;
    }
    void copy_force_point_2(vector <int>& witch)
    {
      for(int ni=0;ni<3;ni++)
	force_point[ni]=(point_ud[witch[0]]->force_point[ni]+point_ud[witch[1]]->force_point[ni])*0.5;
    }
    void copy_force_point_1()
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
    void copy_potential_point_6()
    {
      potential_point=(point_ud[0]->potential_point+point_ud[1]->potential_point+point_ud[2]->potential_point+
		       point_ud[3]->potential_point+point_ud[4]->potential_point+point_ud[5]->potential_point)/6.0;
    }
    void copy_potential_point_4(vector <int>& witch)
    {
      potential_point=(point_ud[witch[0]]->potential_point+point_ud[witch[1]]->potential_point+
		       point_ud[witch[2]]->potential_point+point_ud[witch[3]]->potential_point)*0.25;
    }
    void copy_potential_point_2(vector <int>& witch)
    {
      potential_point=(point_ud[witch[0]]->potential_point+point_ud[witch[1]]->potential_point)*0.5;
    }
    void copy_potential_point_1()
    {
      if(point_pointer != 0)
	{
	  potential_point=point_pointer->potential_point;
	  return;
	}
      dump();
      assert(point_pointer);
    }
    void diff_pot(const double& conv)
    {
      force_point[0]=(get_point_ud_0(0)->potential_point-get_point_ud_0(1)->potential_point)*conv;
      force_point[1]=(get_point_ud_0(2)->potential_point-get_point_ud_0(3)->potential_point)*conv;
      force_point[2]=(get_point_ud_0(4)->potential_point-get_point_ud_0(5)->potential_point)*conv;
    }
    void diff_force(const double& conv)
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
    void div_force(const double& conv)
    {
      // (fx,x)(fx,y)(fx,z)(fy,y)(fy,z)(fz,z)
      force_shear_point.resize(1);
      force_shear_point[0]=
	(-(get_point_ud_0(0)->force_point[0]-get_point_ud_0(1)->force_point[0])
	 -(get_point_ud_0(2)->force_point[0]-get_point_ud_0(3)->force_point[0])
	 -(get_point_ud_0(4)->force_point[0]-get_point_ud_0(5)->force_point[0]))*conv;
    }
    double get_div_force()
    {
      return force_shear_point[0];
    }
    void get_shear(vector <double>& s)
    {
      s=force_shear_point;
    }
    void add_dens_at_point(const double& d)
    {
      density_point+=d;
    }
    template <class T> void add_density_at_points(vector <T>& d)
    {
      if(inside) density_point+=d[0];
      Point* p=get_point_ud_0(1);
      if(p->inside) p->density_point+=d[1];
      p=p->get_point_ud_0(3);
      if(p->inside) p->density_point+=d[3];
      p=p->get_point_ud_0(0);
      if(p->inside) p->density_point+=d[2];
      p=p->get_point_ud_0(5);
      if(p->inside) p->density_point+=d[6];
      p=p->get_point_ud_0(1);
      if(p->inside) p->density_point+=d[7];
      p=p->get_point_ud_0(2);
      if(p->inside) p->density_point+=d[5];
      p=p->get_point_ud_0(0);
      if(p->inside) p->density_point+=d[4];
    }
    void subtract_dens_at_point(const double& d)
    {
      density_point-=d;
    }
    void scale_density_point(const double& s)
    {
      density_point*=s;
    }
    void dens_from_mother()
    {
      density_point=point_pointer->density_point;
    }
    void set_density_point(const double& dens)
    {
      density_point=dens;
    }
    void get_field_values(vector <double>& pott,vector <double>& fx,vector <double>& fy,vector <double>& fz)
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
    void get_field_shear_values(vector <double>& fxx,vector <double>& fxy,vector <double>& fxz,
				vector <double>& fyy,vector <double>& fyz,vector <double>& fzz)
    {
      //      cout << "valuesa " << force_shear_point.size() << endl;
      fxx[0]=force_shear_point[0];
      fxy[0]=force_shear_point[1];
      fxz[0]=force_shear_point[2];
      fyy[0]=force_shear_point[3];
      fyz[0]=force_shear_point[4];
      fzz[0]=force_shear_point[5];
      Point* p=get_point_ud_0(1);
      //      cout << "valuesb " << p->force_shear_point.size() << endl;
      fxx[1]=p->force_shear_point[0];
      fxy[1]=p->force_shear_point[1];
      fxz[1]=p->force_shear_point[2];
      fyy[1]=p->force_shear_point[3];
      fyz[1]=p->force_shear_point[4];
      fzz[1]=p->force_shear_point[5];
      p=p->get_point_ud_0(3);
      //      cout << "valuesc " << p->force_shear_point.size() << endl;
      fxx[3]=p->force_shear_point[0];
      fxy[3]=p->force_shear_point[1];
      fxz[3]=p->force_shear_point[2];
      fyy[3]=p->force_shear_point[3];
      fyz[3]=p->force_shear_point[4];
      fzz[3]=p->force_shear_point[5];
      p=p->get_point_ud_0(0);
      //      cout << "valuesd " << p->force_shear_point.size() << endl;
      fxx[2]=p->force_shear_point[0];
      fxy[2]=p->force_shear_point[1];
      fxz[2]=p->force_shear_point[2];
      fyy[2]=p->force_shear_point[3];
      fyz[2]=p->force_shear_point[4];
      fzz[2]=p->force_shear_point[5];
      p=p->get_point_ud_0(5);
      //      cout << "valuese " << p->force_shear_point.size() << endl;
      fxx[6]=p->force_shear_point[0];
      fxy[6]=p->force_shear_point[1];
      fxz[6]=p->force_shear_point[2];
      fyy[6]=p->force_shear_point[3];
      fyz[6]=p->force_shear_point[4];
      fzz[6]=p->force_shear_point[5];
      p=p->get_point_ud_0(1);
      //      cout << "valuesf " << p->force_shear_point.size() << endl;
      fxx[7]=p->force_shear_point[0];
      fxy[7]=p->force_shear_point[1];
      fxz[7]=p->force_shear_point[2];
      fyy[7]=p->force_shear_point[3];
      fyz[7]=p->force_shear_point[4];
      fzz[7]=p->force_shear_point[5];
      p=p->get_point_ud_0(2);
      //      cout << "valuesg " << p->force_shear_point.size() << endl;
      fxx[5]=p->force_shear_point[0];
      fxy[5]=p->force_shear_point[1];
      fxz[5]=p->force_shear_point[2];
      fyy[5]=p->force_shear_point[3];
      fyz[5]=p->force_shear_point[4];
      fzz[5]=p->force_shear_point[5];
      p=p->get_point_ud_0(0);
      //      cout << "valuesh " << p->force_shear_point.size() << endl;
      fxx[4]=p->force_shear_point[0];
      fxy[4]=p->force_shear_point[1];
      fxz[4]=p->force_shear_point[2];
      fyy[4]=p->force_shear_point[3];
      fyz[4]=p->force_shear_point[4];
      fzz[4]=p->force_shear_point[5];
    }
    void get_field_values(vector <double>& pott)
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
    void get_density_points(vector <double>& dens)
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
    void dumpy()
    {
      cout << "dumpy " << this << " ";
      cout << real_pointer << " ";
      cout << pos_point[0] << " ";
      cout << pos_point[1] << " ";
      cout << pos_point[2] << endl;
    }
    void dumpp()
    {
      cout << "dumpp " << this << " ";
      cout << real_pointer << " " << endl;
      cout << pos_point[0] << " ";
      cout << pos_point[1] << " ";
      cout << pos_point[2] << " ";
      cout << potential_point << endl;
    }
    void dump()
    {
      cout << this << "\t";
      cout << pos_point[0] << "\t";
      cout << pos_point[1] << "\t";
      cout << pos_point[2] << "\t";
      cout << "dump point " << " ";
      cout << real_pointer << endl;
      //
      cout << point_ud[0] << "\t";
      if(point_ud[0] != 0) 
	cout << point_ud[0]->pos_point[0] << "\t" << point_ud[0]->pos_point[1] << "\t" << point_ud[0]->pos_point[2];
      cout << endl;
      cout << point_ud[1] << "\t";
      if(point_ud[1] != 0) 
	cout << point_ud[1]->pos_point[0] << "\t" << point_ud[1]->pos_point[1] << "\t" << point_ud[1]->pos_point[2];
      cout << endl;
      cout << point_ud[2] << "\t";
      if(point_ud[2] != 0) 
	cout << point_ud[2]->pos_point[0] << "\t" << point_ud[2]->pos_point[1] << "\t" << point_ud[2]->pos_point[2];
      cout << endl;
      cout << point_ud[3] << "\t";
      if(point_ud[3] != 0) 
	cout << point_ud[3]->pos_point[0] << "\t" << point_ud[3]->pos_point[1] << "\t" << point_ud[3]->pos_point[2];
      cout << endl;
      cout << point_ud[4] << "\t";
      if(point_ud[4] != 0) 
	cout << point_ud[4]->pos_point[0] << "\t" << point_ud[4]->pos_point[1] << "\t" << point_ud[4]->pos_point[2];
      cout << endl;
      cout << point_ud[5] << "\t";
      if(point_ud[5] != 0) 
	cout << point_ud[5]->pos_point[0] << "\t" << point_ud[5]->pos_point[1] << "\t" << point_ud[5]->pos_point[2];
      cout << endl;
    } 
    void get_field_values(double& pot,double& fx,double& fy,double& fz)
    {
      pot=potential_point;
      fx=force_point[0];
      fy=force_point[1];
      fz=force_point[2];
    }
    void scale_pot_forces(const double& scaling)
    {
      potential_point*=scaling;
      force_point[0]*=scaling;
      force_point[1]*=scaling;
      force_point[2]*=scaling;
    }
    void get_deltas(vector <double>& pos,double& d_x, double& d_y, double& d_z,const double& scale,const double& d_inv)
    {
      //      cout << "deltax " << pos[0] << " " << pos_point[0] << " " << scale << " " << d_inv << endl;
      //      cout << "deltay " << pos[1] << " " << pos_point[1] << " " << scale << " " << d_inv << endl;
      //      cout << "deltaz " << pos[2] << " " << pos_point[2] << " " << scale << " " << d_inv << endl;
      //      cout << " this " << this << endl;
      d_x=(pos[0]*scale-pos_point[0])*d_inv;
      d_y=(pos[1]*scale-pos_point[1])*d_inv;
      d_z=(pos[2]*scale-pos_point[2])*d_inv;
      //      cout << "dx " << d_x << " " << d_y << " " << d_z <<endl;
      assert(abs(d_x-0.5) <= 0.5);
      assert(abs(d_y-0.5) <= 0.5);
      assert(abs(d_z-0.5) <= 0.5);
    }
    Point* neighbor_high(const int& what)
    {
      if(point_ud[what] != 0 && point_ud[what]->get_it_is_high()) return point_ud[what];
      return 0;
    }
    Point* neighbor_high(const int& what1,const int& what2)
    {
      Point* ref=neighbor_high(what1);
      if(ref != 0) return ref->neighbor_high(what2);
      ref=neighbor_high(what2);
      if(ref != 0) return ref->neighbor_high(what1);
      return 0;
    }
    Point* neighbor_high(const int& what1,const int& what2,const int& what3)
    {
      Point* ref=neighbor_high(what1,what2);
      if(ref != 0) return ref->neighbor_high(what3);
      ref=neighbor_high(what2,what3);
      if(ref != 0) return ref->neighbor_high(what1);
      ref=neighbor_high(what1,what3);
      if(ref != 0) return ref->neighbor_high(what2);
      return 0;
    }
    Point* neighbor(const int& what)
    {
      return point_ud[what];
    }
    Point* neighbor(const int& what1,const int& what2)
    {
      Point* ref=neighbor(what1);
      if(ref != 0) return ref->neighbor(what2);
      ref=neighbor(what2);
      if(ref != 0) return ref->neighbor(what1);
      return 0;
    }
    Point* neighbor(const int& what1,const int& what2,const int& what3)
    {
      Point* ref=neighbor(what1,what2);
      if(ref != 0) return ref->neighbor(what3);
      ref=neighbor(what2,what3);
      if(ref != 0) return ref->neighbor(what1);
      ref=neighbor(what1,what3);
      if(ref != 0) return ref->neighbor(what2);
      return 0;
    }
    double sum_inside_weights(vector <double>& w)
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
    inline static bool equal_pos(Point* p1,Point* p2)
    {
      return (p1->pos_point[2]-p2->pos_point[2] == 0) &&
	(p1->pos_point[1]-p2->pos_point[1] == 0) &&
	(p1->pos_point[0]-p2->pos_point[0] == 0);
    }
    inline static bool pointer_test(Point* p)
    {
      return p != 0;
    }
  };
}
#endif
