/*
***  Mr. Garibaldi to hapless Engineer: 
***  It is a prototype, it does not have to be perfect,
***  It just has to work.
     
***  Sheridan to Mr. Garibaldi:
***  I am not interested in your problems, 
***  I am only interested in your solutions.
     
***  Sirian Cybernetics Corporation
***  Its fundamental design flaws are obscured
***  by its superficial design flaws
     
***  This is a prototype of fractal gravity, it is not supposed to be pretty
***  or fast, it is just supposed to work.
***  This is a multi  node version in development
***  with a simplified version of the single node code as a base.
     
***  It "should" be simple to transform into a multi-node code.
***  Code is called repeatedly, it does not remember anything
***  from previous calls except a few useful arrays and the FT of the Green's
***  function for an isolated simulation
     
***  JV Villumsen
***  SHOULD died many years ago.

***   TREESTART
***  generate group 1 at level 0.
***  Now comes the tree generation
***  For each group, find the high points, and arrange them into high groups.
***  Each high group generates one daughter group. Continue this process
***  recursively until we reach level=level_max or 
***  there are no more high points.
***   HIGH_POINTS
***  find all the high points in a group
***   BUFFER_POINTS
***  if necessary, and if we so choose, add buffer high points
***   HIGH_PAIRS
***  find pairs of high points
***   EQUIVALENCE_CLASS
***  find equivalence classes of high points
***   HIGH_GROUPS
***  make high groups from equivalence classes
***   DAUGHTER_GROUP_NINA
***  for each high group generate a daughter group
***  this finishes the tree of groups
***   HIGHEST_LEVEL_GROUP
***  for each particle in the simulation, 
***  find the groups it belongs to and choose the 
***  one at the highest level
***  now loop over all the levels and the groups at the level
***   ASSIGN DENSITY
***  find the density at each point in the group
***   PERIODIC_SOLVER
***   ISOLATED_SOLVER
***  Periodic or isolated FFT solver for potential for group 1 points
***   POTENTIAL_START
***  use the potential values from the mother group to assign initial
***  values to the potential at the group points
***   POISSON_SOLVER
***  solve for the potential in the group 
***   FORCE_AT_POINT
***  difference the potential to get the force field at points
***   FORCE_AT_PARTICLE
***  interpolate to get the forces at the particles that belong to the group
***  CLEAN_UP
***  delete all dynamic memory
***  Now all is done

***  THAT IS ALL FOR NOW FOLKS
*/
#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  int Group::number_groups=0;
  int Point::number_points=0;
  int Particle::number_particles=0;
  bool Point::calc_candidates=true;
  Point* Point::nothing=0;
  bool Fractal::first_time_solver=true;
  void fractal_gravity(Fractal& fractal,Fractal_Memory& fractal_memory)
  {
    //--------------------------------------------------------------------------
    // Generating a class Misc that holds miscellaneous variables and methods invisible to the user
    // Deleted at the end of fractal_gravity
    //--------------------------------------------------------------------------
    if(Point::calc_candidates)
      candidate_points();
    Misc* p_misc=new (nothrow) Misc(fractal.get_tweaks());
    assert(p_misc);
    Misc& misc=*p_misc;
    fractal_memory.p_misc=p_misc;
    cout << "startup " << fractal_memory.start_up << endl;
    misc.zoom=Misc::pow(2,fractal.get_level_max());
    int spam_test=INT_MAX/misc.zoom;
    assert(spam_test > fractal.get_grid_length());
    misc.grid_multiply=misc.zoom*fractal.get_grid_length();  
    misc.set_debug(fractal_memory.debug);
    if(fractal.get_level_start() == 0)
      {
	//--------------------------------------------------------------------------
	// Generating the head group at level 0. It holds a cube of points
	//--------------------------------------------------------------------------
	Group* p_group=new (nothrow) Group;
	assert(p_group);
	Group& group=*p_group;
	group.set_level(0);
	misc.p_group_0=p_group; 
	cout << "group f " << p_group << " " << misc.p_group_0 << endl;
	//--------------------------------------------------------------------------
	// Finding the mean density if periodic BC
	//--------------------------------------------------------------------------
	//	double d_0=0.0; 
	//	if(fractal.get_periodic())
	//	  d_0=total_mass(fractal);
	//	fractal.set_density_0(d_0);
	//--------------------------------------------------------------------------
	fractal_memory.all_groups.resize(fractal.get_level_max()+1);
	fractal_memory.all_groups[0].push_back(p_group);
	//--------------------------------------------------------------------------
	// If isolated BC and this is the first time through fractal_gravity
	// call isolated_solver to generate the FT of the Green's function only.
	//--------------------------------------------------------------------------
	if(Fractal::first_time_solver && !fractal.get_periodic())
	  {
	    fractal.timing(-1,0);
	    isolated_solver(group,fractal,misc); 
	    fractal.timing(1,0);
	  }
	//--------------------------------------------------------------------------
	// In the head group: If periodic BC wrap particles into the unit cube.
	// If isolated BC ignore particles outside unit cube.
	// Generate the points in the head group and assign particles to points.
	//--------------------------------------------------------------------------
	fractal.timing_lev(-2,0);
	fractal.timing(-1,1);
	tree_start(group,fractal,fractal_memory,misc);   
	fractal.timing(1,1);
	fractal.timing_lev(2,0);
      }
    else
      {
	receive_groups(fractal,fractal_memory,misc);   
      }
    //--------------------------------------------------------------------------
    // Generate tree of groups, if depth of tree greater than zero
    //--------------------------------------------------------------------------
    for(int level=fractal.get_level_start();level < fractal.get_level_max();level++)
      {
	fractal.timing_lev(-2,level+1);
	for(vector <Group*>::const_iterator group_itr=fractal_memory.all_groups[level].begin();
	    group_itr!=fractal_memory.all_groups[level].end();group_itr++)
	  {
	    Group& group=**group_itr;
	    if(fractal.get_level_start() == 0 || level > fractal.get_level_start())
	      {
		assert(group.get_number_high_groups() <= 0);
		if(group.get_number_high_groups() == 0) continue;
		//--------------------------------------------------------------------------
		// A high_points is a point with at least "minimum_number" particles (of any mass)
		// associated. The spatial masks specify the highest level at which a point can be a high_point
		//--------------------------------------------------------------------------
		fractal.timing(-1,13);
		high_points(group,fractal,misc);
		fractal.timing(1,13);
		//--------------------------------------------------------------------------
		// If there are any high_points, go on
		//--------------------------------------------------------------------------
		if(group.get_number_high_points() == 0) continue;
		//--------------------------------------------------------------------------
		// If padding or buffering allowed, padding has priority, add buffer points
		// Makes groups bigger, reduces number of spurious small groups, makes boundaries smoother,
		// Ensures that a group never shares a boundary with a daughter group.
		// If padding, all 26 nearest neighbors points of a high_point are made high_points, if not already
		// If buffering, if any of the of the eight octants of the cube defined by the high_points has at least
		// "minimum_number" particles then the seven points nearest to the appropriate corner are made high_points
		//--------------------------------------------------------------------------
		fractal.timing(-1,10);
		buffer_points(group,fractal,misc);
		if(misc.get_debug()) cout << "out of buff_points" << endl;
		group.p_list_really_high.clear();
		fractal.timing(1,10);
		//--------------------------------------------------------------------------
		// Two points are a pair if, and only if, the two cubes they define have a face in common
		// A point generates a max of three pairs.
		//--------------------------------------------------------------------------
		fractal.timing(-1,11);
		if(misc.get_debug()) cout << "into high_pairs" << endl;
		high_pairs(group,misc); 
		if(misc.get_debug()) cout << "out of high_pairs" << endl;
		fractal.timing(1,11);
		//--------------------------------------------------------------------------
		// Find equivalence classes of high_points
		//--------------------------------------------------------------------------
		fractal.timing(-1,12);
		if(misc.get_debug()) cout << "into eq class" << endl;
		equivalence_class(group);
		if(misc.get_debug()) cout << "out of eq class" << endl;
		group.list_pair_1.clear();
		group.list_pair_2.clear();
		fractal.timing(1,12);
		//--------------------------------------------------------------------------
		// Groups of high_points make high_groups as equivalence classes of high_points
		//--------------------------------------------------------------------------
		fractal.timing(-1,13);
		if(misc.get_debug()) cout << "into high_groups" << endl;
		high_groups(group);
		if(misc.get_debug()) cout << "out of high_groups" << endl;
		split_high_groups(group,fractal,fractal_memory,misc);
		group.head_number.resize(0);
		group.list_high.resize(0);  	
		fractal.timing(1,13);
		send_groups(group,fractal,fractal_memory,misc);
	      }
	    //--------------------------------------------------------------------------
	    // Loop over all high_groups
	    //--------------------------------------------------------------------------
	    for(vector <Group*>::const_iterator high_group_itr=group.list_high_groups.begin();
		high_group_itr != group.list_high_groups.end();++high_group_itr)
	      {
		Group& high_group=**high_group_itr;
		//--------------------------------------------------------------------------
		// Each all high_group makes one new group one level higher
		//--------------------------------------------------------------------------
		Group* p_new_group=new (nothrow) Group(group);
		assert(p_new_group);
		Group& new_group=*p_new_group;
		new_group.set_mother_group(&group);
		//--------------------------------------------------------------------------
		// Generate the points in the new group, get rid of duplicates, find the six nearest
		// neighbor points, if they exist. Decide if a point is "inside", assign particles to points
		//--------------------------------------------------------------------------
		if(misc.get_debug()) cout << "into daughter" << endl;
		daughter_group_nina(new_group,high_group,fractal,fractal_memory,misc);
		if(misc.get_debug()) cout << "out of daughter" << endl;
		//--------------------------------------------------------------------------
		// High_group no longer needed, no memory leaks please
		// Never instantiate "point.p_in_high_group", it is only a label.
		//--------------------------------------------------------------------------
		delete &high_group;
		//--------------------------------------------------------------------------
		// Decide if the new group can have daughter groups
		//--------------------------------------------------------------------------
		if(new_group.get_points_in_group() == 27)
		  new_group.set_number_high_groups(0);
		else
		  new_group.set_number_high_groups(-1);
		//--------------------------------------------------------------------------
		fractal_memory.all_groups[level+1].push_back(p_new_group);
	      }
	  }
	fractal.timing_lev(2,level+1);
      }
    fractal.timing(-1,16);
    bool badd=test_tree(fractal_memory,fractal);
    assert(!badd);
    fractal.timing(1,16);
    //--------------------------------------------------------------------------
    // If force_max > 0 find ghost particles that need to split to soften the force
    //--------------------------------------------------------------------------
    fractal.timing(-1,17);
    Fractal* p_fractal_ghost=new (nothrow) Fractal;
    assert(p_fractal_ghost);
    Fractal& fractal_ghost=*p_fractal_ghost;
    heavies(fractal,fractal_ghost,misc);
    fractal.timing(1,17);
    //--------------------------------------------------------------------------
    // Each ghost particle, if any, is assigned to a point in each group it belongs to
    //--------------------------------------------------------------------------
    fractal.timing(-1,18);
    particle_lists(fractal_memory.all_groups,fractal,fractal_ghost,misc);
    fractal.timing(1,18);
    //--------------------------------------------------------------------------
    if(!fractal_memory.start_up)
      {
	fractal.timing_lev(-1,0);
	for(int level=fractal.get_level_start();level <= fractal.get_level_max();level++)
	  {
	    fractal.timing_lev(-1,level);
	    for(vector <Group*>::const_iterator group_itr=fractal_memory.all_groups[level].begin();
		group_itr!=fractal_memory.all_groups[level].end();group_itr++)
	      {
		Group& group=**group_itr;
		if(misc.get_debug())
		  cout << "trying doing Group " << &group << " " << group.get_level() <<  endl;
		group.set_set_scaling(true);
		group.set_set_dens(true);
		//--------------------------------------------------------------------------
		// use all particles to assign density to the points using cic. For edge points, density
		//  will be wrong, but this density is never used.
		//--------------------------------------------------------------------------
		if(level == 0)
		  {
		    fractal.timing(-1,3);
		    assign_density(group,fractal); 
		    fractal.timing(1,3);
		    if(fractal.get_periodic())
		      {
			//--------------------------------------------------------------------------
			// For periodic BC use periodic FFTW solver for head group
			//--------------------------------------------------------------------------
			fractal.timing(-1,4);
			periodic_solver(group,fractal_memory,fractal); 
			fractal.timing(1,4);
		      }
		    else
		      {
			//--------------------------------------------------------------------------
			// For isolated BC use isolated FFTW solver for head group
			//--------------------------------------------------------------------------
			fractal.timing(-1,5);
			isolated_solver(group,fractal,misc); 
			fractal.timing(1,5);
		      }
		  }
		else
		  {
		    fractal.timing(-1,19);
		    assign_density(group,fractal); 
		    fractal.timing(1,19);
		    //--------------------------------------------------------------------------
		    // For all other groups use mother group potential as initial
		    // conditions and boundary conditions for potential
		    //--------------------------------------------------------------------------
		    receive_potential(group,fractal,fractal_memory,misc);
		    fractal.timing(-1,20);
		    potential_start(group); 		
		    fractal.timing(1,20);
		    //--------------------------------------------------------------------------
		    // Use SOR or some other Poisson solver for the potential
		    //--------------------------------------------------------------------------
		    fractal.timing(-1,21);
		    poisson_solver(group,fractal,misc); 		
		    fractal.timing(1,21);
		    send_potential(group,fractal,fractal_memory,misc);
		  }
		//--------------------------------------------------------------------------
		// For inside points diff potential to get forces at points. For all other points
		// use values from mother group
		//--------------------------------------------------------------------------
		fractal.timing(-1,22);
		force_at_point(group,fractal); 		   
		fractal.timing(1,22);
		//--------------------------------------------------------------------------
		// If this group is the highest level group for this particle use cic interpolation
		// to find potential and forces at the particle
		//--------------------------------------------------------------------------
		fractal.timing(-1,23);
		force_at_particle(group,fractal);
		group.set_done_group(true);
		fractal.timing(1,23);
	      }
	    fractal.timing_lev(1,level);
	  }
      }
    if(fractal.get_halo_fixed()) 
      {
	fractal.timing(-1,25);
	halo_force_fixed(fractal);
	fractal.timing(1,25);
      }
    if(fractal_memory.start_up)
      {
	cout << " going to particle_initial " << endl;
	initial_forces_sharp(fractal_memory,fractal,misc);
      }
    if(fractal_memory.calc_shear || fractal_memory.calc_density_particle)
      {
	for(int level=fractal.get_level_start();level <= fractal.get_level_max();level++)
	  {
	    cout << "level = " << level << endl;
	    for(vector <Group*>::const_iterator group_itr=fractal_memory.all_groups[level].begin();
		group_itr!=fractal_memory.all_groups[level].end();group_itr++)
	      {
		force_shear_at_point(**group_itr,fractal);
	      }
	  }
	cout << "enter shear at part " << endl;
	force_shear_at_particle(fractal_memory,fractal,misc);
	cout << "exit shear at part " << endl;
      }
    //--------------------------------------------------------------------------
    // Clean up all dynamically allocated memory except fractal
    //--------------------------------------------------------------------------
    groups_level(fractal,fractal_memory.all_groups);
    fractal.timing(-1,26);
    clean_up(fractal_memory,misc,fractal_ghost);
    fractal.timing(1,26);
  }
}
