/*
***  Mr. Garibaldi to hapless Engineer: 
***  It is a prototype, it does not have to be perfect,
***   it just has to work.
     
***  Sheridan to Mr. Garibaldi:
***  I am not interested in your problems, 
***   I am only interested in your solutions.
     
***  Sirian Cybernetics Corporation:
***  Share and Enjoy.
***  Its fundamental design flaws are obscured
***   by its superficial design flaws.

***  JV Villumsen:
***  The Answer to Life, the Universe and Everything is 42.
***   My dorm room number at Aarhus University was 42.
***    To quote Douglas Adams: "I thinks that's terribly significant".
***  Should died many years ago.
***  If all else fails read the instructions.
***  !the_Messiah = a_naughty_boy.
***  If it were easy, it would have been done long time ago.
***  You do not have to be faster than the bear. You just have
***   to be faster than the other guy running from the bear.
***  If I knew what I was doing, it would not be science,
***   it would be engineering.

***  Sheldon Cooper:
***  I am a Physicist,
***   I have a working Knowledge of the Universe and Everything in it.

***  Bomb #20: Intriguing. I wish I had more time to discuss this.
***  Doolittle: [frantic] Why don't you have more time?
***  Bomb #20: Because I must explode in 75 seconds. 

***  Bomb, this is Lt. Doolittle. You are *not* to detonate in the bomb bay. I repeat, you are NOT to detonate in the bomb bay! 

***  KISS
***  Keep It Simple Stupid

***  This is a prototype of fractal force, it is not supposed to be pretty
***  or fast, it is just supposed to work.
***  This is a multi node version
     
***  Code is called repeatedly, it does not remember anything
***  from previous calls except a few useful arrays and the FT of the Green's
***  function for an isolated simulation

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
***  now loop over all the chains and the groups in a chain
***   ASSIGN DENSITY
***  find the density at each point in the group
***   PERIODIC_SOLVER
***   ISOLATED_SOLVER
***  Periodic or isolated FFT solver for potential for group 0 points
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
  void fractal_force(Fractal& fractal,Fractal_Memory& fractal_memory)
  {
    cout << "here in fractal_force " << endl;
    scatter_particles(fractal_memory,fractal);

    if(Point::calc_candidates)
      candidate_points();
    if(fractal_memory.mind_wipe)
      {
	Fractal* p_fractal_dummy=new (nothrow) Fractal;
	assert(p_fractal_dummy);
	Fractal& fractal_dummy=*p_fractal_dummy;
	clean_up(fractal_memory,*(fractal_memory.p_misc),fractal_dummy);
	fractal_memory.p_misc=0;
	fractal_memory.all_groups.clear();
	return;
      }
    if(fractal_memory.fixed_potential)
      {
	fractal.timing(-1,18);
	particle_lists_fixed(fractal_memory.all_groups,fractal,*(fractal_memory.p_misc));
	fractal.timing(1,18);
	fractal.timing(-1,23);
	force_at_particle(fractal_memory.all_groups,fractal);
	fractal.timing(1,23);
	if(fractal.get_halo_fixed()) 
	  {
	    fractal.timing(-1,25);
	    halo_force_fixed(fractal);
	    fractal.timing(1,25);
	  }
	return;
      }
    Misc* p_misc=new (nothrow) Misc;
    assert(p_misc);
    Misc& misc=*p_misc;
    fractal_memory.p_misc=p_misc;
    cout << "startup " << fractal_memory.start_up << endl;
    misc.zoom=Misc::pow(2,fractal.get_level_max());
    misc.grid_multiply=misc.zoom*fractal.get_grid_length();  
    misc.set_debug(fractal_memory.debug);
    //--------------------------------------------------------------------------
    // Finding the mean density if periodic BC
    //--------------------------------------------------------------------------
    double d_0=0.0; 
    if(fractal.get_periodic())
      d_0=fractal.get_omega_start()*0.375/(4.0*atan(1.0));
    fractal.set_density_0(d_0);
    //--------------------------------------------------------------------------
    // Generating the head group at level 0. It holds a cube of points
    //--------------------------------------------------------------------------
    Group* p_group=new (nothrow) Group;
    assert(p_group);
    Group& group=*p_group;
    group.set_level(0);
    misc.p_group_0=p_group; 
    cout << "group f " << p_group << " " << misc.p_group_0 << endl;
    fractal_memory.all_groups.resize(fractal.get_level_max()+1);
    fractal_memory.all_groups[0].push_back(p_group);
    //--------------------------------------------------------------------------
    // If isolated BC and this is the first time through fractal_force
    // call isolated_solver to generate the FT of the Green's function only.
    //--------------------------------------------------------------------------
    if(Fractal::first_time_solver && !fractal.get_periodic())
      {
	fractal.timing(-1,0);
	isolated_solver(group,fractal_memory,fractal,misc); 
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
    //
    //
    //--------------------------------------------------------------------------
    // Generate tree of groups, if depth of tree greater than zero
    //--------------------------------------------------------------------------
    // Loop over all levels to generate new groups at the next level
    // if there are high_points in the group. Done recursively until
    // no group generates high_points. The n'th recursion generates groups at level "n".
    // A group does not generate a new group if (1) no high_points, (2) group is at level=level_max
    // (3) Group only has 27 points
    //--------------------------------------------------------------------------
    for(int level=0;level < fractal.get_level_max();level++)
      {
	fractal.timing_lev(-2,level+1);
	for(vector <Group*>::const_iterator group_itr=fractal_memory.all_groups[level].begin();
	    group_itr!=fractal_memory.all_groups[level].end();group_itr++)
	  {
	    Group& group=**group_itr;
	    assert(group.get_number_high_groups() <= 0);
	    if(group.get_number_high_groups() == 0) continue;
	    //--------------------------------------------------------------------------
	    // A high_points is a point with at least "minimum_number" particles (of any mass)
	    // associated. The spatial masks specify the highest level at which a point can be a high_point
	    //--------------------------------------------------------------------------
	    fractal.timing(-1,9);
	    high_points(group,fractal,misc);
	    fractal.timing(1,9);
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
	    fractal.timing(1,10);
	  }
	//
	for(vector <Group*>::const_iterator group_itr=fractal_memory.all_groups[level].begin();
	    group_itr!=fractal_memory.all_groups[level].end();group_itr++)
	  {
	    Group& group=**group_itr;
	    assert(group.get_number_high_groups() <= 0);
	    if(group.get_number_high_groups() == 0) continue;
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
	    fractal.timing(1,12);
	    //--------------------------------------------------------------------------
	    // Groups of high_points make high_groups as equivalence classes of high_points
	    //--------------------------------------------------------------------------
	    fractal.timing(-1,13);
	    if(misc.get_debug()) cout << "into high_groups" << endl;
	    high_groups(group);
	    if(misc.get_debug()) cout << "out of high_groups" << endl;
	    group.head_number.resize(0);
	    group.list_high.resize(0);  	
	    fractal.timing(1,13);
	    //--------------------------------------------------------------------------
	    // Loop over all high_groups
	    //--------------------------------------------------------------------------
	    int group_counter=0;
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
		new_group.set_group_number(group_counter);
		//--------------------------------------------------------------------------
		// Generate the points in the new group, get rid of duplicates, find the six nearest
		// neighbor points, if they exist. Decide if a point is "inside", assign particles to points
		//--------------------------------------------------------------------------
		//		if(misc.get_debug()) cout << "into daughter" << endl;
		daughter_group_nina(new_group,high_group,fractal,fractal_memory,misc);
		//		if(misc.get_debug()) cout << "out of daughter" << endl;
		//--------------------------------------------------------------------------
		// High_group no longer needed, no memory leaks please
		// Never instantiate "point.p_in_high_group", it is only a label.
		//--------------------------------------------------------------------------
		delete &high_group;
		//--------------------------------------------------------------------------
		// Decide if the new group can have daughter groups
		//--------------------------------------------------------------------------
		if(new_group.get_points_in_group() == 27 && !fractal_memory.MPIrun)
		  new_group.set_number_high_groups(0);
		else
		  new_group.set_number_high_groups(-1);
		//--------------------------------------------------------------------------
		group_counter++;
		fractal_memory.all_groups[level+1].push_back(p_new_group);
	      }
	  }
	fractal.timing_lev(2,level+1);
      }
    //    dump_tree(fractal_memory,fractal);
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
    if(!fractal_memory.start_up)
      {
	// solve level zero by first assigning density the use either periodic or isolated FFTW
	// to find potential at the points. Find forces at points and then at particles.
	fractal.timing_lev(-1,0);
	Group& group=*misc.p_group_0;
	if(misc.get_debug())
	  cout << "trying doing Group " << &group << " " << group.get_level() <<  endl;
	group.set_set_scaling(true);
	group.set_set_dens(true);
	//--------------------------------------------------------------------------
	// use all particles to assign density to the points using cic. For edge points, density
	//  will be wrong, but this density is never used. Also add the density from any split ghost particles
	//--------------------------------------------------------------------------
	fractal.timing(-1,3);
	//			if(fractal_memory.calc_density_particle && group.get_level() > 0) density_edge(group,misc);
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
	    isolated_solver(group,fractal_memory,fractal,misc); 
	    fractal.timing(1,5);
	  }
	//--------------------------------------------------------------------------
	// For inside points diff potential to get forces at points. For all other points
	// use values from mother group
	//--------------------------------------------------------------------------
	fractal.timing(-1,7);
	force_at_point(group,fractal); 		   
	fractal.timing(1,7);
	//--------------------------------------------------------------------------
	// If this group is the highest level group for this particle use cic interpolation
	// to find potential and forces at the particle
	//--------------------------------------------------------------------------
	fractal.timing(-1,8);
	force_at_particle(group,fractal,fractal_memory.momentum_conserve);
	group.set_done_group(true);
	fractal.timing(1,8);
	fractal.timing_lev(1,0);
	// loop over all levels > 0 
	for(int level=1;level <= fractal.get_level_max();level++)
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
		//  will be wrong, but this density is never used. Also add the density from any split ghost particles
		//--------------------------------------------------------------------------
		fractal.timing(-1,19);
		assign_density(group,fractal); 
		fractal.timing(1,19);
		//--------------------------------------------------------------------------
		// For all other groups use mother group potential as initial
		// conditions and boundary conditions for potential
		//--------------------------------------------------------------------------
		fractal.timing(-1,20);
		potential_start(group); 		
		fractal.timing(1,20);
	      }
	    //--------------------------------------------------------------------------
	    // Use SOR or some other Poisson solver for the potential
	    //--------------------------------------------------------------------------
	    fractal.timing(-1,21);
	    poisson_solver(fractal,fractal_memory,level); 		
	    fractal.timing(1,21);
	    for(vector <Group*>::const_iterator group_itr=fractal_memory.all_groups[level].begin();
		group_itr!=fractal_memory.all_groups[level].end();group_itr++)
	      {
		Group& group=**group_itr;
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
		if(fractal_memory.momentum_conserve)
		  force_sum_particle(group,fractal,false);
		fractal.timing(-1,23);
		force_at_particle(group,fractal,fractal_memory.momentum_conserve);
		if(fractal_memory.momentum_conserve)
		  force_sum_particle(group,fractal,true);
		group.set_done_group(true);
		fractal.timing(1,23);
	      }
	    fractal.timing_lev(1,level);
	  }
      }
    if(fractal_memory.momentum_conserve)
      {
	force_test(fractal_memory,fractal);
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
	for(int level=0;level <= fractal.get_level_max();level++)
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
    gather_particles(fractal_memory,fractal);
  }
}
