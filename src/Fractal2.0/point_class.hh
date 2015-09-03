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
    bool found_it;
    bool really_passive;
    bool passive_point;
    bool buffer_point;
    bool edge_point;
    bool mass_point;
    int number_in_list;
    //    int which_Slice;
    //
    int ij_number;
    vector <int>ij_ud;
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
    //    vector <Particle*> list_other_particles;
    static ofstream* p_FILE;
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
      found_it(false),
      really_passive(false),
      passive_point(false),
      buffer_point(false),
      edge_point(false),
      mass_point(false),
      number_in_list(-1),
      potential_point(0.0),
      density_point(0.0),
      pos_point(3,-1),
      point_ud(6,nothing),
      force_point(3,0.0)
    {
      number_points++;
    }
    ~Point()
    {    
      number_points--;
    }
    void set_found_it(bool fi);
    bool get_found_it() const;
    void set_mass_point(bool what);
    void set_mass_points(bool what);
    bool get_mass_point() const;
    void set_edge_buffer_passive_point(const bool& e,const bool& b,const bool& p);
    void set_edge_buffer_passive_really_point(const bool& e,const bool& b,const bool& p,const bool& r);
    bool get_buffer_point() const;
    void set_buffer_point(bool& b_p);
    bool get_edge_point() const;
    void set_edge_point(bool& e_p);
    int get_number_in_list() const;
    void set_number_in_list(const int& n);
    void set_eureka_adj(const vector <bool>& eu);
    void get_eureka_adj(vector <bool>& eu) const;
    void set_eureka_dau(const vector <bool>& eu);
    void get_eureka_dau(vector <bool>& eu) const;
    bool get_eureka_dau(const int& ni) const;
    void set_p_in_group(Group* p_g);
    Group* get_p_in_group() const;
    void set_p_in_high_group(Group* p_g);
    Group* get_p_in_high_group();
    void get_pos_point(vector <int>& pos) const;
    void set_pos_point(const vector <int>& pos);
    void set_pos_point(const int& x,const int& y,const int& z);
    void get_pos_point(int& x,int& y,int& z) const;
    int get_pos_point_x() const;
    void set_pos_point_x(const int& i);
    int get_pos_point_y() const;
    void set_pos_point_y(const int& i);
    int get_pos_point_z() const;
    void set_pos_point_z(const int& i);
    int get_pos_point(const int& i) const;
    int get_real_pointer() const;
    void set_real_pointer(const int& i);
    void set_point_to_number(const int& i);
    int get_point_to_number() const;
    void set_high_number(const int& i);
    int get_high_number() const;
    bool get_inside() const;
    void set_inside(const bool& ins);
    bool get_passive_point() const;
    void set_passive_point(const bool& value);
    bool get_really_passive() const;
    void set_really_passive(const bool& value);
    void set_ij_number(const int& count);
    int get_ij_number() const;
    void really_clear(vector <Point*>& die);
    void set_ij_neighbors();
    void set_ij_neighbors(vector <int>& Box);
    void get_ij_neighbors(vector <int>& ijud) const;
    int get_ij_neighbors_size() const;
    void copy_ij_index(const int& ijc);
    void get_hypre_info(int& ij_index,vector <int>& ijud,double& rho,double& pot) const;
    bool get_it_is_high() const;
    void set_it_is_high(const bool& value);
    void set_passive_low();
    void set_inside_high();
    double get_potential_point() const;
    void set_potential_point(const double& d);
    void add_potential_point(const double& d);
    void add_force_point_x(const double& d);
    void add_force_point_y(const double& d);
    void add_force_point_z(const double& d);
    double get_density_point() const;
    void get_force_point(vector <double>& force) const;
    void set_force_point(vector <double>& force);
    double get_force_point_x() const;
    double get_force_point_y() const;
    double get_force_point_z() const;
    void set_force_point_x(const double& f);
    void set_force_point_y(const double& f);
    void set_force_point_z(const double& f);
    void set_FILE(ofstream* p_filE);
    ofstream* get_FILE();
//     void all_mine(vector <Point*>& pointers,vector <bool>& belongs_to_me);
    Point* move_adj(const int& ra,const int& rb);
    Point* move_rp(const int& r);
    void get_point_ud(vector <Point*>& point_6) const;
    void set_point_ud(vector <Point*>& point_6);
    Point* get_point_up_x() const;
    Point* get_point_up_x_0() const;
    void set_point_up_x(Point* p_point);
    Point* get_point_up_y() const;
    Point* get_point_up_y_0() const;
    void set_point_up_y(Point* p_point);
    Point* get_point_up_z() const;
    Point* get_point_up_z_0() const;
    void set_point_up_z(Point* p_point);
    Point* get_point_down_x() const;
    Point* get_point_down_x_0() const;
    void set_point_down_x(Point* p_point);
    void down_from_up(Point* p_up_x,Point* p_up_y,Point* p_up_z);
    void down_from_up();
    void up_from_down();
    Point* get_point_down_y() const;
    Point* get_point_down_y_0() const;
    void set_point_down_y(Point* p_point);
    Point* get_point_down_z() const;
    Point* get_point_down_z_0() const;
    void set_point_down_z(Point* p_point);
    Point* get_point_ud(const int& i) const;
    Point* get_point_ud_0(const int& i,const int& tag) const;
    Point* get_point_ud_0(const int& i) const;
    void set_point_ud(Point* p_point,const int& i);
    Point* get_point_pointer() const;
    Point* get_point_pointer_t() const;
    void set_point_pointer(Point* p_point);
    Point* get_p_daughter_point() const;
    void set_p_daughter_point(Point* p_d_point);
    void point_pointers_all(Point& high_point);
    double laplacian() const;
    void force_shear_point_zero();
    void copy_force_shear_point(Point& p0);
    void copy_density_point(Point& p0,Point& p1);
    void copy_density_point(Point& p0,Point& p1,Point& p2,Point& p3);
    void copy_potential_point(Point& p0);
    void copy_force_shear_point_6();
    void copy_force_shear_point_4(vector <int>& witch);
    void copy_force_shear_point_2(vector <int>& witch);
    void copy_force_shear_point_1();
    void copy_force_point_6();
    void copy_force_point_4(vector <int>& witch);
    void copy_force_point_2(vector <int>& witch);
    void copy_force_point_1();
    void copy_potential_point_6();
    void copy_potential_point_4(vector <int>& witch);
    void copy_potential_point_2(vector <int>& witch);
    void copy_potential_point_1();
    void diff_pot(const double& conv);
    void diff_pot_careful(const double& conv);
    void diff_force(const double& conv);
    void div_force(const double& conv);
    double get_div_force() const;
    void get_shear(vector <double>& s) const;
    void add_dens_at_point(const double& d);
    template <class T> void add_density_at_points(vector <T>& d);
    void subtract_dens_at_point(const double& d);
    void scale_density_point(const double& s);
    void dens_from_mother();
    void set_density_point(const double& dens);
    void get_field_values(vector <double>& pott,vector <double>& fx,vector <double>& fy,vector <double>& fz) const;
    void get_field_shear_values(vector <double>& fxx,vector <double>& fxy,vector <double>& fxz,
				       vector <double>& fyy,vector <double>& fyz,vector <double>& fzz) const;
    void get_field_values(vector <double>& pott) const;
    void get_density_points(vector <double>& dens);
    void dumpy() const;
    void dumpd() const;
    void dumpp() const;
    void dumpp(ofstream& FF) const;
    void dumppf() const;
    void dump() const;
    void get_field_values(double& pot,double& fx,double& fy,double& fz) const;
    void scale_pot_forces(const double& scaling);
    void get_deltas(vector <double>& pos,double& d_x, double& d_y, double& d_z,const double& scale,const double& d_inv) const;
    void get_potss_denss(const int* dprl,const double& g_c,double* potss,double* denss);
    void set_potss(const int* dprl,const double* potss);
    Point* neighbor_high(const int& what) const;
    Point* neighbor_high(const int& what1,const int& what2) const;
    Point* neighbor_high(const int& what1,const int& what2,const int& what3) const;
    Point* neighbor(const int& what) const;
    Point* neighbor(const int& what1,const int& what2) const;
    Point* neighbor(const int& what1,const int& what2,const int& what3) const;
    double sum_inside_weights(vector <double>& w);
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
