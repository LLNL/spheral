namespace FractalSpace
{
  void add_buffer_values(Fractal_Memory& mem,int level,
			 vector < vector<int> >& SBoxes,vector < vector<Point*> >& SPoints);
  void add_pseudo_particles(Fractal_Memory& mem,Fractal& frac);
  void adj_nina(Point& point,vector <Point*>& adj);
  double Age_of_the_universe (const double& omega_0, const double& omega_lambda, const double& redshift);
  template <class ForwardIterator>
  void am_I_conservative_enough_isol(Fractal_Memory* PFM,ForwardIterator massesb,double G,
				     vector <double>& xmin,vector <double>& xmax,double correction,
				     ForwardIterator posxb,ForwardIterator posyb,
				     ForwardIterator poszb,ForwardIterator velxb,
				     ForwardIterator velyb,ForwardIterator velzb);
  // void am_I_conservative_enough_isol(Fractal_Memory* PFM,vector <double>& masses,double G,
  // 				     vector <double>& xmin,vector <double>& xmax,double correction,
  // 				     vector <double>& posx,vector <double>& posy,vector <double>& posz,
  // 				     vector <double>& velx,vector <double>& vely,vector <double>& velz);
  void any_overlaps(Fractal_Memory& mem,int spacing,int VOLMIN,double FILLFACTOR,vector<vector<int>>& SBoxes,vector<vector<Point*>>& SPoints);
  void assign_density(Group& group, Fractal& fractal);
  void balance_by_particles(Fractal_Memory* PFM,bool withparts);
  void binary_balancing(vector <double>& numbers,double minimum,
			int Nodes,int length,vector <double>& targets,vector <int>& lowers,vector <int>& uppers);
  void binary_balancing(Fractal_Memory* PFM,vector <double>& numbers,double minimum,
			int Nodes,int length,vector <double>& targets,vector <int>& lowers,vector <int>& uppers);
  void box_stats(Fractal_Memory& mem,int level,int nb,vector<vector<int>>& SBoxes,vector<vector<Point*>>& SPoints);
  void buffer_points(Group& group, Fractal& fractal);
  void candidate_points();
  void check_for_edge_trouble(Fractal& fractal);
  bool check_high(Point& point,Fractal& fractal);
  bool check_high(Point& point,Fractal_Memory& fractal_memory);
  void clean_groups(Fractal_Memory& fractal_memory);
  void clean_shear(Fractal_Memory& fractal_memory);
  void clean_overlaps(Fractal_Memory& mem,int spacing,int VOLMIN,double FILLFACTOR,vector<bool>& STrouble,vector<vector<int>>& SBoxes,vector<vector<Point*>>& SPoints);
  void clean_up(Fractal_Memory& mem,Misc& misc);
  template <typename T> void clean_vector(vector<T>& vec);
  template <typename T> void clean_deque(deque<T>& deq);
  bool compare_vectorsX(vector <int> veca,vector <int> vecb);
  bool compare_vectorsY(vector <int> veca,vector <int> vecb);
  bool compare_vectorsZ(vector <int> veca,vector <int> vecb);
  double cosmos_power(const double& k,Fractal_Memory& fractal_memory);
  void daughter_group_nina(Group& new_group,Group& high_group,Fractal& fractal,Fractal_Memory& memo,Misc& misc);
  void density_edge(Group& group, Misc& misc);
  void dens_to_slices(Group& group,Fractal_Memory& mem,Fractal& frac);
  template <class T> T det(const vector <T>& m);
  double dGrowthdA(const double& omega_0, const double& omega_lambda, const double& redshift);
  double dGrowthdT(const double& omega_0, const double& omega_lambda, const double& redshift);
  double dist1(const double& x,const double& y);
  void dump(Point& point);
  void dump_all_particles(Fractal& fractal);
  void dump_cosmo_boxes(const int& step,Fractal& fractal);
  void dump_group(Group& group,Misc& misc);
  void dump_tree(Fractal_Memory& fractal_memory,Fractal& fractal);
  void edge_buffer_inside(vector <int>& n,vector <int>& Box,vector <int>& BBox,vector <int>& Buffer,bool& MPIrun,
			  vector <bool>& periods,bool& inside,vector <bool>& buff,vector <bool>& edge);
  template <class M>  void energy_simple(M& mem, Fractal& fractal);
  void equivalence_class(Group& group);
  void factors(int FR,vector <int>& divs,bool& easy);
  void find_global_level_max(Fractal_Memory& mem);
  int findPointByPos(vector <Point*>& plist,Point* psend,FILE* PFF);
  void fix_memory(Fractal& frac,const int& ispace,const int& jfield);
  void force_at_particle(Group& group, Fractal& fractal,const bool& conserve);
  void force_at_particle(vector <vector <Group*> >& all_groups, Fractal& fractal);
  void force_at_particle_sharp(Group& group, Fractal& fractal);
  void force_at_point(Group& group, Fractal& fractal);
  void force_at_point_initial(Group& group, Fractal& fractal);
  void force_shear_at_particle(Fractal_Memory& fractal_memory,Fractal& fractal);
  void force_shear_at_point(Group& group,Fractal& fractal);
  void force_sum_particle(Group& group,const bool& doit);
  void force_test(Fractal& fractal);
  void fractal_force(Fractal& fractal,Fractal_Memory& fractal_memory);
  void fractal_force_init(Fractal_Memory* p_mem);
  int fractal_force_wrapper(Fractal_Memory* p_fractal_memory,Fractal* p_fractal);
  template <class M> void fractal_memory_parameters(M* pmem,const string _disK_,int _mulT_);
  void fractal_memory_parameters(Fractal_Memory* pmem,double _mulT_);
  void fractal_nina_cosmo(int argc,char* argv[]);
  void fractal_nina_galaxy(int argc,char* argv[]);
  void Full_Stop(Fractal_Memory& mem,int number);
  void Full_Stop(Fractal_Memory& mem,MPI_Comm& Comm,int number);
  void gather_particles(Fractal_Memory& mem,Fractal& frac);
  void go_ahead_points(vector <Point*>& adj,vector <bool>& ins,vector <bool>& go_ahead);
  bool group_in_box(Group* pgroup,vector<int>& BOX);
  void groups_level(Fractal& fractal,vector < vector<Group*> >& all_groups);
  double Growth(const double& omega_0, const double& omega_lambda, const double& redshift);
  void halo_force(Fractal& fractal);
  void halo_force_fixed(Fractal& frac);
  void heavies(Fractal& fractal,Fractal& fractal_ghost);
  void highest_level_group(Fractal_Memory& fractal_memory,Fractal& fractal,Misc& misc);
  void highest_level_group_other(vector <vector <Group*> >& all_groups,Fractal& fractal,Misc& misc);
  bool high_enough_level(Point& point,Group& group,Fractal& fractal,Misc& misc);
  void high_groups(Group& group);
  void high_pairs(Group& group);
  void high_points(Group& group, Fractal& fractal,Misc& misc);
  double Hubble (const double& omega_0, const double& omega_lambda, const double& redshift);
  void hypre_best_boxes(Fractal_Memory& mem,vector<vector<Point*> >& hypre_points,int spacing,int& VOLbest,double& FILLbest);
  void hypre_best_boxes(bool buffer,Fractal_Memory& mem,vector<vector<Point*> >& hypre_points,int spacing,int& VOLbest,double& FILLbest);
  void hypre_clever_boxes(Fractal_Memory& mem,vector <vector<Point*>>& hypre_points,int spacing,
			  int VOLMIN,double FILLFACTOR,
			  vector < vector<int> >& SBoxes,vector < vector<Point*> >& SPoints);
  void hypre_dump(int level,vector <Point*>& hypre_points,ofstream& FH);
  void hypre_eror(FILE* PFH,int level,int ni,int er);
  bool hypre_struct_load_balance(Fractal_Memory& mem,vector<vector<int>>& SBoxes,vector<vector<Point*>>& SPoints,vector<int>& HRout);
  void hypre_points_boxes(Fractal_Memory& mem,vector <vector <Point*> >& hypre_points,int spacing,
			  int VOLMIN,double FILLFACTOR,
			  vector < vector<int> >& SBoxes,vector < vector<Point*> >& SPoints);
  void hypre_points_clean(Fractal_Memory& mem,int level,vector< vector<Point*> >& hypre_points);
  void hypre_points_struct(bool single,Fractal_Memory& mem,vector <Group*>& groups,
			   vector < vector <Point*> >& hypre_points,bool buffer_groups,int level);
  void hypre_points_zero(vector<vector<Point*>>& SPoints);
  void hypre_solve_struct(bool buffer,Fractal_Memory& mem,int level,
			  vector < vector<int> >& SBoxes,vector < vector<Point*> >& SPoints);
  void hypre_test_boxes(Fractal_Memory& mem,int level,
			vector < vector<int> >& SBoxes,vector < vector<Point*> >& SPoints);
  void hypre_world_create(Fractal_Memory& mem,int level,vector <vector <int> >& SBoxes,vector<vector<Point*>>& SPoints,
			  bool buffer_groups);
  void hypre_world_destroy();
  void info_to_slices(Fractal_Memory& mem,Fractal& frac,int lev);
  void info_to_slices_to_pot_init(Fractal_Memory& mem,Fractal& frac,int lev);
  void initial_forces_sharp(Fractal_Memory& fractal_memory,Fractal& fractal);
  void isolated_solver(Group& group,Fractal_Memory& fractal_memory,Fractal& fractal);
  bool it_is_inside(Point* p_point);
  inline bool it_is_inside(Point* p_point);
  double Lambda (const double& omega_0, const double& omega_lambda, const double& redshift);
  double Length(const double& omega_0, const double& omega_lambda, const double& redshift);
  bool LesserPoint(Point* p1,Point* p2);
  bool LesserPointA(Point* p1,Point* p2);
  void left_right(deque <Point*>& all_points,vector <int>& pos_left,vector <int>& pos_right);
  void left_right(deque <Point*>& all_points,vector <int>& pos_left,vector <int>& pos_right,const int& wrap);
  void left_right(vector <Point*>& all_points,vector <int>& pos_left,vector <int>& pos_right);
  void left_right(vector <Point*>& all_points,vector <int>& pos_left,vector <int>& pos_right,const int& wrap);
  void left_right(Fractal& frac,vector <double>& pos_left,vector <double>& pos_right);
  void left_right(vector <Group*>& all_groups,vector <int>& pos_left,vector <int>& pos_right);
  void list_buffer(Point& point,const int& corner);
  template <class ForwardIterator>
  void make_me_a_galaxy(Fractal_Memory* PFM,int numbers,double total_mass,double G,
			ForwardIterator posxb,ForwardIterator posyb,ForwardIterator poszb,
			ForwardIterator velxb,ForwardIterator velyb,ForwardIterator velzb);
  void make_me_a_galaxy(int FractalRank,int numbers,double total_mass,double G,
			vector <double>& xpos,vector <double>& ypos,vector <double>& zpos,
			vector <double>& xvel,vector <double>& yvel,vector <double>& zvel);
  void make_me_some_particles(int rank,int numbers,double total_mass,vector <double>& masses,double G,
			      vector <double>& xpos,vector <double>& ypos,vector <double>& zpos,
			      vector <double>& xvel,vector <double>& yvel,vector <double>& zvel);
  void make_decisions_erika(Misc& misc);
  template <class M, class F>  void make_particles(M& mem,F& frac,int& count,const double& m,const bool& crash);
  void match_edges(Fractal_Memory& mem,int level);
  void max_predict(Fractal_Memory& fractal_memory,Fractal& fractal,vector <double>& shear_force,double& min_vol);
  bool mini_solve(Fractal_Memory& mem,Group* pg);
  void mini_solve1(Point* p,const double& gc);
  void mini_solve2(Point* pa,Point* pb,const double& gc);
  void mini_solve3(const vector<Point*>& found,const double& gc);
  void move_small_boxes(Fractal_Memory& mem,vector<int>& Boxes,vector<vector<Point*>>& SPoints,vector<int>& HRout);
  void neighbor_easy(vector <Point*>& p);
  void neighbors_nina(Point& point, vector <Point*>& adj);
  void node_groups_struct(Fractal_Memory& mem,vector <int>& counts);
  double Omega (const double& omega_0, const double& omega_lambda, const double& redshift);
  bool on_edge(vector <int>& pos,vector <int>& Box);
  bool on_edge(array <int,3>& pos,vector <int>& Box);
  bool on_edge(Point* p,vector <int>& Box);
  template <class T> bool overlap(vector <T>& xleft,vector <T>& xright,vector <T>& yleft,vector <T>& yright);
  template <class T> bool overlap(vector <T>& xleft,vector <T>& xright,vector <T>& box);
  template <class T> bool overlap_boxes(vector <T>& xvec,vector <T>& box);
  template <class T> bool overlap_interval(T Imin,T Imax,T Jmin,T Jmax,T& LOW,T& HIGH);
  void particle_lists(vector <vector <Group*> >& all_groups,Fractal& fractal,Fractal& fractal_ghost,Misc& misc);
  void particle_lists_fixed(vector <vector <Group*> >& all_groups,Fractal& fractal,Misc& misc);
  void periodic_solver(Group& group, Fractal_Memory& fractal_memory,Fractal& fractal);
  void points_on_nodes(Fractal_Memory& mem);
  void poisson_solver_struct(Fractal& fractal,Fractal_Memory& mem,const int& level);
  void potential_start(Group& group);
  void power_spectrum(fftw_complex* rhoC,int length,vector <double>& variance_rho,vector <double>& variance_pot,
		      vector <double>& variance_force,vector <double>& variance_force_s,int lev,double d0,bool do_var,
		      Fractal_Memory& mem);
  bool rad_compare(Particle* par1,Particle* par2);
  bool right_diff(vector <int>& Va,vector <int>& Vb,vector <int>& VD);
  void remove_dupe_points(int spacing,vector<vector<Point*>>& hypre_points,vector<vector<int>>& SBoxes,vector<vector<Point*>>& SPoints);
  void remove_pseudo_particles(Fractal_Memory& mem,Fractal& frac);
  void scatter_particles(Fractal_Memory& mem,Fractal& frac);
  template <class T> int shortest_vector(vector<T>& veca,vector<T>& vecb,vector<T>& vecc);
  template <class ForwardIterator>
  void shrink_cube(double SHRINK,
		   const vector <double>& xmin,const vector <double>& xmax,
		   Fractal_Memory* PFM,
		   ForwardIterator posxb,
		   ForwardIterator posyb,
		   ForwardIterator poszb,
		   int number_particles,
		   vector <double>& xmini,vector <double>& xmaxy);
  void shrink_cube(Fractal_Memory* PFM,double SHRINK,const vector <double>& xmin,const vector <double>& xmax,
		   vector <double>& xmini,vector <double>& xmaxy);
  void shrink_cube(double SHRINK,const vector <double>& xmin,const vector <double>& xmax,Fractal_Memory* PFM,
		   vector <double>& posx,vector <double>& posy,vector <double>& posz,
		   int number_particles,vector <double>& xmini,vector <double>& xmaxy);
  void slices_to_potf(Fractal_Memory& mem,Fractal& frac,int lev);
  void slices_to_pot_init(Fractal_Memory& mem,Fractal& frac,int lev);
  void small_exceptions(Fractal_Memory& mem);
  void sor(Group& group, Fractal& fractal,vector <Point*>& list_left_x,const int& dir);
  void sor_solver(Group& group, Fractal& fractal);
  void sort3_list(Group& group,int what);
  void sort3_list(vector <Point*>& list_points,int what);
  void sort3_list(deque <Point*>& list_points,int what);
  void sort_3(Fractal& fractal,Group& group);
  int spawn(const vector <double>& pos, vector < vector <double> > ppos,const vector <double>& boxouter);
  void split_nodes(int FR,int& FR0,int& FR1,int& FR2);
  template <class M, class F> int split_particle(M& mem,F& frac,const double& x0,const double& y0,const double& z0,
						 int& count,const double& m,const int& split_to,const bool& gen_part);
  template <class ForwardIterator>
  void start_writing(Fractal_Memory* PFM,int Numberparticles,
		     double G,vector <double>& xmin,vector <double>& xmax,
		     ForwardIterator posxb,ForwardIterator posyb,ForwardIterator poszb,
		     ForwardIterator velxb,ForwardIterator velyb,ForwardIterator velzb,ForwardIterator massesb);
  void start_writing(Fractal_Memory* PFM,int Numberparticles,double G,vector <double>& xmin,vector <double>& xmax,
		     vector<double>& posx,vector<double>& posy,vector<double>& posz,
		     vector<double>& velx,vector<double>& vely,vector<double>& velz,vector<double>& masses);
  template <class M>  void step_simple(M& mem,Fractal& fractal);
  void sum_pot_forces(Fractal& fractal);
  void super_groups(Fractal_Memory& mem,vector <Group*>& groups,const int level,
		    vector<vector<int>>& WorldRanks,
		    vector<vector<int>>& LocalGroups,
		    vector<vector<int>>& FreeNodes,
		    vector<bool>& IAmIn);
  template <class ForwardIterator>
  void take_a_leap_isol(Fractal_Memory* PFM,ForwardIterator massesb,double G,
			vector <double>& xmin,vector <double>& xmax,
			ForwardIterator posxb,ForwardIterator posyb,
			ForwardIterator poszb,ForwardIterator velxb,
			ForwardIterator velyb,ForwardIterator velzb);
  void take_a_leap_isol(Fractal_Memory* PFM,vector <double>& masses,double G,
			vector <double>& xmin,vector <double>& xmax,
			vector <double>& posx,vector <double>& posy,vector <double>& posz,
			vector <double>& velx,vector <double>& vely,vector <double>& velz);
  void test_gal(Fractal_Memory& mem,Fractal& fractal);
  bool test_good_point(Point* p1,Fractal_Memory& mem,int level);
  bool test_group(Group& group);
  void test_points(Fractal_Memory& mem,vector<vector<Point*>>& SPoints,int level);
  bool test_tree(Fractal_Memory& fractal_memory,Fractal& fractal);
  void tree_dump(Fractal_Memory& FM);
  void tree_start(Group& group,Fractal& fractal,Fractal_Memory& memo,Misc& misc);
  void tree_start_mini(Group& group,Fractal& fractal,Fractal_Memory& memo,Misc& misc);
  Point* try_harder(Point& point0,const int& ni,const bool& easy);
  void update_rv(Fractal& fractal,const int& param,const double& const1,const double& const2);
  void use_freenodes(Fractal_Memory& mem,vector<int>& countsP,vector<int>& countsB);
  template <class T> bool vector_in_box(vector <T>& xvec,vector <T>& box);
  template <class T> bool vector_in_box(array <T,3>& xvec,vector <T>& box);
  bool vector_in_box(Point* p,vector <int>& box);
  bool vector_in_box(const Point& p,vector <int>& box);
  void velocities(Fractal_Memory& mem,Fractal& frac);
  int which_element(vector <Point*>& vec,int x,int y,int z,bool periodic,int period,ofstream& FF);
  void write_rv(const int& step,Fractal& fractal);
}
