/*
***  Mr. Garibaldi to hapless Engineer: 
***  It is a prototype, it does not have to be perfect,
***   it just has to work.
     
***  Sheridan to Mr. Garibaldi:
***  I am not interested in your problems, 
***   I am only interested in your solutions.
     
***  Sirius Cybernetics Corporation:
***   Share and Enjoy.
***  Its fundamental design flaws are obscured
***   by its superficial design flaws.

***  Urban Dictionary:
***   Stupiphany
***   A sudden, intuitive perception of or insight into the reality or essential meaning of something very basic, 
***   annoyingly obvious, or - in hindsight - really stupid. 
***   The person having the stupiphany is generally greatly excited, while observers just shake their heads at how obvious the answer was.

***  JV Villumsen:
***  The Answer to Life, the Universe and Everything is 42.
***   My dorm room number at Aarhus University was 42.
***    Douglas Adams: "I thinks that's terribly significant".
***  Should died many years ago.
***  If all else fails read the instructions.
***  !the_Messiah = a_naughty_boy.
***  If it were easy, it would have been done long time ago.
***  You do not have to be faster than the bear. You just have
***   to be faster than the other guy running from the bear.
***  If I knew what I was doing, it would not be science,
***   it would be engineering.
***  In college I brought the computer system for Western Denmark to it's knees,
***   with a FORTRAN program.
***  For my PhD thesis I used 1200 particles on a single node. PDP 11/60.
***   Now I run 1000000 particles per node on 100000 nodes. Vulcan@LLNL.
***   Eight orders of magnitude.
***  On an observing trip in Chile I was interrogated at gun-point, in Spanish, about my suspected involvement
***   in an armed bank robbery.
***  On my honeymoon our delicious meal was interrupted by a shooting in the restaurant.
***   The owner limped over to our table and explained that he had been shot by the cook.
***   He apologized for the incident, said he was going to the hospital, but asked us to please enjoy our meal.
***  I need to find out where I can get a Pan Galactic Gargle Blaster.
***  TC Villumsen at 13 y.o.:
***   Mom, Dad does not do any work, all he does is read magazines, play on the computer and talk with his friends.

***  Mike Owen:
***  Our simulations are better than reality.
***  As soon as you get a program to work 
***   you should rewrite it from scratch.

***  Sheldon Cooper:
***  I am a Physicist,
***   I have a working Knowledge of the Universe and Everything in it.

***  Bomb #20: Intriguing. I wish I had more time to discuss this.
***  Doolittle: [frantic] Why don't you have more time?
***  Bomb #20: Because I must explode in 75 seconds. 

***  Bomb, this is Lt. Doolittle. You are *not* to detonate in the bomb bay. I repeat, you are NOT to detonate in the bomb bay! 

***  KISS
***  Keep It Simple Stupid

*** The Big Bang Theory
*** Our whole universe was in a hot dense state,
*** Then nearly fourteen billion years ago expansion started. Wait...
*** The Earth began to cool,
*** The autotrophs began to drool,
*** Neanderthals developed tools,
*** We built a wall (we built the pyramids),
*** Math, science, history, unraveling the mystery,
*** That all started with the big bang!
*** Since the "Dawn of Man" is really not that long,
*** As every galaxy was formed in less time than it takes to sing this song.
*** A fraction of a second and the elements were made.
*** The bipeds stood up straight,
*** The dinosaurs all met their fate,
*** They tried to leave but they were late
*** And they all died (they froze their asses off)
*** The oceans said, "Pangaea,
*** See ya, wouldn't wanna be ya!"
*** Set in motion by the same big bang!
*** It all started with the big BANG!
*** It's expanding ever outward but one day
*** It will pause then start to go the other way,
*** Collapsing ever inward, we won't be here, it won't be heard.
*** Our best and brightest figure that it'll make an even bigger bang!
*** Australopithecus would really have been sick of us,
*** Debating how we're here; they're catching deer (we're catching viruses)
*** Religion or astronomy, Descartes or Deuteronomy
*** It all started with the big bang!
*** Music and mythology, Einstein and astrology
*** It all started with the big bang!
*** It all started with the big ... BANG! 

*** Curiosity Rover
***  Tango Delta Nominal
***  RIMU Stable
***  UHF Strong

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
***  now loop over all the levels and the groups at that level
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
    Full_Stop(fractal_memory,-1);
    ofstream& FileFractal=fractal.p_file->FileFractal;
    FileFractal << "here in fractal_force " << endl;
    FileFractal << "number of everything entering fractal "  << " " << Group::number_groups << " " << Point::number_points << endl;
    FileFractal << " Total number of particles entering Fractal " << Particle::number_particles << endl;
    //    write_rv(-16,fractal);
    Full_Stop(fractal_memory,33);
    fractal.timing_lev(-1,0);
    fractal.timing(-1,27);
    scatter_particles(fractal_memory,fractal);
    fractal.timing(1,27);
    fractal.timing_lev(1,0);
    //    write_rv(-15,fractal);
    int jfield=4;
    if(fractal_memory.calc_density_particle) 
      jfield=5;
    fix_memory(fractal,3,jfield);
    if(Point::calc_candidates)
      candidate_points();
    Misc* p_misc=new (nothrow) Misc;
    assert(p_misc);
    Misc& misc=*p_misc;
    bool debug=fractal.get_debug();
    misc.set_debug(debug);
    fractal_memory.p_misc=p_misc;
    FileFractal << "startup " << fractal_memory.start_up << endl;
    misc.zoom=Misc::pow(2,fractal.get_level_max());
    misc.grid_multiply=misc.zoom*fractal.get_grid_length();  
    misc.set_debug(fractal_memory.debug);
    //--------------------------------------------------------------------------
    // Finding the mean density if periodic BC
    //--------------------------------------------------------------------------
    double d_0=0.0; 
    if(fractal.get_periodic())
      {
	d_0=fractal.get_omega_start()*0.375/(4.0*atan(1.0));
	FileFractal << "density 0 " << fractal.get_omega_start() << " " << d_0 << endl;
      }
    fractal.set_density_0(d_0);
    //--------------------------------------------------------------------------
    // Generating the head group at level 0. It holds a box of points
    //--------------------------------------------------------------------------
    Group* p_group=new (nothrow) Group;
    assert(p_group);
    Group& group=*p_group;
    group.set_level(0);
    group.p_file=fractal.p_file;
    misc.p_group_0=p_group; 
    FileFractal << "group f " << p_group << " " << misc.p_group_0 << endl;
    fractal_memory.all_groups.resize(fractal.get_level_max()+1);
    fractal_memory.all_groups[0].push_back(p_group);
    fractal_memory.all_buffer_groups.resize(fractal.get_level_max()+1);
    fractal_memory.all_buffer_groups[0].push_back(p_group);
    fractal_memory.all_inside_groups.resize(fractal.get_level_max()+1);
    //--------------------------------------------------------------------------
    // If isolated BC and this is the first time through fractal_force
    // call isolated_solver to generate the FT of the Green's function only.
    //--------------------------------------------------------------------------
    if(Fractal::first_time_solver && !fractal.get_periodic())
      {
	fractal.timing(-1,0);
	isolated_solver(group,fractal_memory,fractal); 
	fractal.timing(1,0);
      }
    //--------------------------------------------------------------------------
    // In the head group: If periodic BC wrap particles into the unit cube.
    // If isolated BC ignore particles outside unit cube.
    // Generate the points in the head group and assign particles to points.
    //--------------------------------------------------------------------------
    fractal.timing(-1,46);
    fractal.timing_lev(-2,0);
    fractal.timing(-1,1);
    tree_start(group,fractal,fractal_memory,misc);   
    fractal.timing(1,1);
    fractal.timing_lev(2,0);
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
	FileFractal << "a level=  " << level << endl;
	fractal.timing_lev(-2,level+1);
	for(vector <Group*>::const_iterator group_itr=fractal_memory.all_groups[level].begin();
	    group_itr!=fractal_memory.all_groups[level].end();group_itr++)
	  {
	    Group& group=**group_itr;
	    //	    FileFractal << "a group=  " << *group_itr << endl;
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
	    //	    FileFractal << "aa group=  " << *group_itr << " " << group.get_number_high_points() << endl;
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
	    if(misc.get_debug()) FileFractal << "out of buff_points" << endl;
	    fractal.timing(1,10);
	  }
	//
	int group_counter=0;
	for(vector <Group*>::const_iterator group_itr=fractal_memory.all_groups[level].begin();
	    group_itr!=fractal_memory.all_groups[level].end();group_itr++)
	  {
	    Group& group=**group_itr;
	    FileFractal << "b group=  " << *group_itr << " " << group.get_number_high_groups() << endl;
	    assert(group.get_number_high_groups() <= 0);
	    if(group.get_number_high_groups() == 0) continue;
	    //--------------------------------------------------------------------------
	    // Two points are a pair if, and only if, the two cubes they define have a face in common
	    // A point generates a max of three pairs.
	    //--------------------------------------------------------------------------
	    fractal.timing(-1,11);
	    FileFractal << "into high_pairs" << endl;
	    high_pairs(group); 
	    FileFractal << "out of high_pairs" << endl;
	    fractal.timing(1,11);
	    //--------------------------------------------------------------------------
	    // Find equivalence classes of high_points
	    //--------------------------------------------------------------------------
	    fractal.timing(-1,12);
	    if(misc.get_debug()) FileFractal << "into eq class" << endl;
	    equivalence_class(group);
	    if(misc.get_debug()) FileFractal << "out of eq class" << endl;
	    fractal.timing(1,12);
	    //--------------------------------------------------------------------------
	    // Groups of high_points make high_groups as equivalence classes of high_points
	    //--------------------------------------------------------------------------
	    fractal.timing(-1,13);
	    if(misc.get_debug()) FileFractal << "into high_groups" << endl;
	    high_groups(group);
	    if(misc.get_debug()) FileFractal << "out of high_groups" << group.list_high_groups.size() << endl;
	    group.head_number.resize(0);
	    group.list_high.resize(0);  	
	    fractal.timing(1,13);
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
		new_group.set_group_number(group_counter);
		//--------------------------------------------------------------------------
		// Generate the points in the new group, get rid of duplicates, find the six nearest
		// neighbor points, if they exist. Decide if a point is "inside", assign particles to points
		//--------------------------------------------------------------------------
		if(misc.get_debug()) FileFractal << "into daughter" << endl;
		daughter_group_nina(new_group,high_group,fractal,fractal_memory,misc);
		if(misc.get_debug()) FileFractal << "out of daughter" << endl;
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
		if(p_new_group->get_buffer_group())
		  fractal_memory.all_buffer_groups[level+1].push_back(p_new_group);
		else
		  fractal_memory.all_inside_groups[level+1].push_back(p_new_group);
	      }
	  }
	fractal.timing_lev(2,level+1);
      }
    //    dump_tree(fractal_memory,fractal);
    fractal.timing(-1,16);
    bool badd=test_tree(fractal_memory,fractal);
    FileFractal << " not badd " << badd << endl;
    assert(!badd);
    fractal.timing(1,16);
    //--------------------------------------------------------------------------
    // If force_max > 0 find ghost particles that need to split to soften the force
    //--------------------------------------------------------------------------
    fractal.timing(-1,17);
    FileFractal << " not badd0 " << badd << endl;
    Fractal* p_fractal_ghost=new (nothrow) Fractal;
    assert(p_fractal_ghost);
    Fractal& fractal_ghost=*p_fractal_ghost;
    FileFractal << " not badda " << badd << endl;
    heavies(fractal,fractal_ghost);
    FileFractal << " not baddb " << badd << endl;
    fractal.timing(1,17);
    //--------------------------------------------------------------------------
    // Each ghost particle, if any, is assigned to a point in each group it belongs to
    //--------------------------------------------------------------------------
    fractal.timing(-1,18);
    FileFractal << " not baddc " << badd << endl;
    particle_lists(fractal_memory.all_groups,fractal,fractal_ghost,misc);
    FileFractal << " not baddd " << badd << endl;
    fractal.timing(1,18);
    fractal.timing(1,46);
    //    fractal_memory.p_mess->TreeTime=fractal.get_delta_time(46);
    FileFractal << "start up " << fractal_memory.start_up << endl;
    if(!fractal_memory.start_up)
      {
	fractal.timing(-1,47);
	// solve level zero by first assigning density then use either periodic or isolated FFTW
	// to find potential at the points. Find forces at points and then at particles.
	fractal.timing_lev(-1,0);
	Group& group=*misc.p_group_0;
	if(misc.get_debug())
	  FileFractal << "trying doing Group " << &group << " " << group.get_level() <<  endl;
	group.set_set_scaling(true);
	group.set_set_dens(true);
	//--------------------------------------------------------------------------
	// use all particles to assign density to the points using cic. For edge points, density
	//  will be wrong, but this density is never used. Also add the density from any split ghost particles
	//--------------------------------------------------------------------------
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
	    isolated_solver(group,fractal_memory,fractal); 
	    fractal.timing(1,5);
	  }
	Full_Stop(fractal_memory,41);
	fractal_memory.global_level_max=find_global_level_max(fractal_memory,fractal);
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
	for(int level=1;level <= fractal_memory.global_level_max;level++)
	  {
	    fractal.timing_lev(-1,level);
	    for(vector <Group*>::const_iterator group_itr=fractal_memory.all_groups[level].begin();
		group_itr!=fractal_memory.all_groups[level].end();group_itr++)
	      {
		Group& group=**group_itr;
		if(misc.get_debug())
		  FileFractal << "trying doing Group " << &group << " " << group.get_level() <<  endl;
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
	    Full_Stop(fractal_memory,36);
	    fractal.timing(-1,31);
	    poisson_solver(fractal,fractal_memory,level); 		
	    fractal.timing(1,31);
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
		  force_sum_particle(group,false);
		fractal.timing(-1,23);
		force_at_particle(group,fractal,fractal_memory.momentum_conserve);
		if(fractal_memory.momentum_conserve)
		  force_sum_particle(group,true);
		group.set_done_group(true);
		fractal.timing(1,23);
	      }
	    fractal.timing_lev(1,level);
	  }
	fractal.timing(1,47);
      }
    if(fractal_memory.momentum_conserve)
      {
	FileFractal << " going to mom conserve " << endl;
	force_test(fractal);
	FileFractal << " leaving mom conserve " << endl;
      }
    if(fractal_memory.start_up)
      {
	fractal.timing(-1,45);
	FileFractal << " going to particle_initial " << endl;
	initial_forces_sharp(fractal_memory,fractal);
	fractal.timing(1,45);
      }
    if(fractal_memory.calc_shear || fractal_memory.calc_density_particle)
      {
	for(int level=0;level <= fractal_memory.global_level_max;level++)
	  {
	    FileFractal << "level = " << level << endl;
	    for(vector <Group*>::const_iterator group_itr=fractal_memory.all_groups[level].begin();
		group_itr!=fractal_memory.all_groups[level].end();group_itr++)
	      {
		force_shear_at_point(**group_itr,fractal);
	      }
	  }
	FileFractal << "enter shear at part " << endl;
	force_shear_at_particle(fractal_memory,fractal);
	FileFractal << "exit shear at part " << endl;
      }
    //--------------------------------------------------------------------------
    // Clean up all dynamically allocated memory except fractal
    //--------------------------------------------------------------------------
    groups_level(fractal,fractal_memory.all_groups);
    if(fractal_memory.steps % fractal_memory.number_steps_out == 0)
      tree_dump(fractal_memory);
    fractal.timing(-1,26);
    clean_up(fractal_memory,misc,fractal_ghost);
    fractal.timing(1,26);
    //    write_rv(-14,fractal);
    Full_Stop(fractal_memory,38);
    fractal.timing_lev(-1,0);
    fractal.timing(-1,28);
    gather_particles(fractal_memory,fractal);
    fractal.timing(1,28);
    fractal.timing_lev(1,0);
    //    write_rv(-13,fractal);
    if(fractal.get_halo_fixed()) 
      {
	fractal.timing(-1,25);
	halo_force_fixed(fractal);
	fractal.timing(1,25);
      }
    //    if(!fractal.get_periodic() && debug)
    //          test_gal(fractal_memory,fractal);
    //    assert(0);
    Full_Stop(fractal_memory,-1);
    FileFractal << "number of everything exiting Fractal "  << " " << Group::number_groups << " " << Point::number_points << endl;
    FileFractal << " Total number of particles exiting Fractal " << Particle::number_particles << endl;
    FileFractal << " Made It fractal_force " << fractal_memory.steps << " " << fractal_memory.p_mess->Clock()-fractal_memory.p_mess->WallTime << endl;
  }
  void Full_Stop(Fractal_Memory& mem,int number)
    {
      Fractal* pf=mem.p_fractal;
      if(number >= 0)
      	pf->timing(-1,number);
      mem.p_mess->Full_Stop();
      if(number >= 0)
      	pf->timing(1,number);
    }
}
