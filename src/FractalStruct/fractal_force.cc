/*
***  Rasmus Villumsen
***   In Honor to the forgotten polar explorer.

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
***  !the_Messiah == a_naughty_boy.
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
***  On an observing trip in Chile I was interrogated at gun-point, in Spanish, 
***   about my suspected involvement in an armed bank robbery.
***  In Brazil I took a picture of the Danish flag on the building housing the Danish consulate.
***   This building also housed Chase Manhattan offices, so I was detained by security for an hour. Shakedown.
***  On my honeymoon our delicious meal was interrupted by a shooting in the restaurant.
***   The owner limped over to our table and explained that he had been shot by the cook.
***   He apologized for the incident, said he was going to the hospital, but asked us to please enjoy our meal.
***  I need to find out where I can get a Pan Galactic Gargle Blaster.
***  In customs at JFK they became suspicious of me, so I was given an astronomy quiz. I passed.
***  To my elementary school teacher, Mrs Hansen, a great teacher. Quiet authority.
***  To my 7th grade teacher, Mr Christensen, a former sailor, who told the best stories, some of them G rated. He made me want to travel the world.
***  To my high school history teacher, Mr. Frederiksen, who taught me critical thinking. Always question everything.
***  To my high school math teacher, Mr. Mourier, who taught me problem solving skills. Infinite patience.
***  TC Villumsen at 13 y.o.:
***   Mom, Dad does not do any work, all he does is read magazines, play on the computer and talk with his friends.

***  Nina:
***  I am not a kid, I am a highly sophisticated child.
***  GrandPa, you are weird.

***  Mike Owen:
***  Our simulations are better than reality..
***  As soon as you get a program to work 
***   you should rewrite it from scratch.

***  Big Bang Theory: Sheldon Cooper:
***  I am a Physicist,
***   I have a working Knowledge of the Universe and Everything in it.

***  Big Bang Theory: Stuart:
***  Like shooting nerds in a barrel

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

***  Bomb #20: Intriguing. I wish I had more time to discuss this.
***  Doolittle: [frantic] Why don't you have more time?
***  Bomb #20: Because I must explode in 75 seconds. 

***  Bomb, this is Lt. Doolittle. You are *not* to detonate in the bomb bay. I repeat, you are NOT to detonate in the bomb bay! 

***  KISS
***  Keep It Simple Stupid

***  TopGun
***  I feel the need, the urgent need for speed

*** Curiosity Rover
***  Tango Delta Nominal
***  RIMU Stable
***  UHF Strong

***  Toast(On Hawaii Five-0)
***  Your nerd credentials have been verified.

***  This is a prototype of fractal force, it is not supposed to be pretty
***  or fast, it is just supposed to work.
     
***  Code is called repeatedly, it does not remember anything
***  from previous calls except a few useful arrays and the FT of the Green's
***  function for an isolated simulation

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
    ofstream& FileFractal=fractal.p_file->DUMPS;
    FileFractal << "here in fractal_force " << "\n";
    FileFractal << "number of everything entering fractal "  << " " << Group::number_groups << " " << Point::number_points << "\n";
    FileFractal << " Total number of particles entering Fractal " << Particle::number_particles << "\n";
    Full_Stop(fractal_memory,33);
    fractal.timing_lev(-1,0);
    scatter_particles(fractal_memory,fractal);
    vector <int> PBoxLength(3);
    fractal.getPBoxLength(PBoxLength);
    int volume=PBoxLength[0]*PBoxLength[1]*PBoxLength[2];
    FileFractal << " particles after scatter " << fractal_memory.p_mess->FractalRank << " " << fractal.get_number_particles() << " " << volume << "\n";
    fractal.timing_lev(1,0);
    int jfield=4;
    if(fractal_memory.calc_density_particle) 
      jfield=5;
    fix_memory(fractal,3,jfield);
    if(Point::calc_candidates)
      candidate_points();
    Misc* p_misc=new Misc;
    assert(p_misc);
    Misc& misc=*p_misc;
    bool debug=fractal.get_debug();
    misc.set_debug(debug);
    fractal_memory.p_misc=p_misc;
    misc.zoom=Misc::pow(2,fractal.get_level_max());
    misc.grid_multiply=misc.zoom*fractal.get_grid_length();  
    misc.set_debug(fractal_memory.debug);
    double d_0=0.0; 
    if(fractal.get_periodic())
      d_0=fractal.get_omega_start()*0.375/(4.0*atan(1.0));
    fractal.set_density_0(d_0);
    Group* p_group=new Group;
    assert(p_group);
    Group& group=*p_group;
    group.set_level(0);
    group.p_file=fractal.p_file;
    misc.p_group_0=p_group; 
    fractal_memory.all_groups.resize(fractal.get_level_max()+1);
    fractal_memory.all_groups[0].push_back(p_group);
    Point::p_FILE=&(fractal_memory.p_mess->p_file->DUMPS);
    fractal.timing(-1,46);
    fractal.timing_lev(-2,0);
    if(fractal.get_periodic())
      tree_start(group,fractal,fractal_memory,misc);   
    else
      tree_start_mini(group,fractal,fractal_memory,misc);   
    fractal.timing_lev(2,0);
    if(!fractal_memory.start_up)
      {
	fractal.timing(-1,47);
	fractal.timing_lev(-1,0);
	Group& group=*misc.p_group_0;
	group.set_set_scaling(true);
	group.set_set_dens(true);
	/*                */	         fractal.timing(-1,3);
	assign_density(group,fractal); 
	/*                */	         fractal.timing(1,3);
	if(fractal.get_periodic())
	  periodic_solver(group,fractal_memory,fractal); 
	else
	  isolated_solver(group,fractal_memory,fractal); 
	/*                */      	fractal.timing(-1,7);
	force_at_point(group,fractal); 		   
	/*                */    	fractal.timing(1,7);
	/*                */    	fractal.timing(-1,8);
	force_at_particle(group,fractal,fractal_memory.momentum_conserve);
	/*                */    	fractal.timing(1,8);
	group.set_done_group(true);
	fractal.timing_lev(1,0);
      }
    for(int level=0;level < fractal.get_level_max();level++)
      {
	fractal.timing_lev(-2,level+1);
	fractal_memory.level=level;
	for(Group* pgroup : fractal_memory.all_groups[level])
	  {
	    Group& group=*pgroup;
	    assert(group.get_number_high_groups() <= 0);
	    if(group.get_number_high_groups() == 0) continue;
	    high_points(group,fractal,misc);
	    if(group.get_number_high_points() == 0) continue;
	    buffer_points(group,fractal);
	  }
	int group_counter=0;
	for(Group* pgroup : fractal_memory.all_groups[level])
	  {
	    Group& group=*pgroup;
	    assert(group.get_number_high_groups() <= 0);
	    if(group.get_number_high_groups() == 0) continue;
	    fractal.timing(-1,11);
	    high_pairs(group); 
	    fractal.timing(1,11);
	    fractal.timing(-1,12);
	    equivalence_class(group);
	    fractal.timing(1,12);
	    fractal.timing(-1,13);
	    high_groups(group);
	    fractal.timing(1,13);
	    for(Group* phigh_group : group.list_high_groups)
	      {
		Group& high_group=*phigh_group;
		Group* p_new_group=new Group(group);
		assert(p_new_group);
		Group& new_group=*p_new_group;
		new_group.set_mother_group(&group);
		new_group.set_group_number(group_counter);
		daughter_group_nina(new_group,high_group,fractal,fractal_memory,misc);
		delete &high_group;
		if(new_group.get_points_in_group() == 27 && !fractal_memory.MPIrun)
		  new_group.set_number_high_groups(0);
		else
		  new_group.set_number_high_groups(-1);
		group_counter++;
		fractal_memory.all_groups[level+1].push_back(p_new_group);
	      }
	  }
	fractal.timing_lev(2,level+1);
      }
    bool badd=test_tree(fractal_memory,fractal);
    assert(!badd);
    FileFractal << "number of everything after the tree "  << " " << Group::number_groups << " " << Point::number_points << "\n";
    FileFractal << " Total number of particles after the tree " << Particle::number_particles << "\n";
    // Fractal* p_fractal_ghost=new Fractal;
    // Fractal& fractal_ghost=*p_fractal_ghost;
    fractal.timing(1,46);
    if(!fractal_memory.start_up)
      {
	Full_Stop(fractal_memory,41);
	find_global_level_max(fractal_memory);
	points_on_nodes(fractal_memory);
	for(int level=1;level <= fractal_memory.global_level_max;level++)
	  {
	    fractal.timing_lev(-1,level);
	    fractal_memory.level=level;
	    for(Group* pgroup : fractal_memory.all_groups[level])
	      {
		Group& group=*pgroup;
		group.set_set_scaling(true);
		group.set_set_dens(true);
	/*                */	         fractal.timing(-1,19);
		assign_density(group,fractal); 
	/*                */	         fractal.timing(1,19);
	/*                */	         fractal.timing(-1,20);
		potential_start(group); 		
	/*                */	         fractal.timing(1,20);
	      }
	    Full_Stop(fractal_memory,36);
	    poisson_solver_struct(fractal,fractal_memory,level);
	    for(Group* pgroup : fractal_memory.all_groups[level])
	      {
		Group& group=*pgroup;
		fractal.timing(-1,22);
		force_at_point(group,fractal); 		   
		fractal.timing(1,22);
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
      force_test(fractal);
    if(fractal_memory.start_up)
      {
	fractal.timing(-1,45);
	FileFractal << " going to particle_initial " << "\n";
	initial_forces_sharp(fractal_memory,fractal);
	fractal.timing(1,45);
      }
    if(fractal_memory.calc_shear || fractal_memory.calc_density_particle)
      {
	for(int level=0;level <= fractal_memory.global_level_max;level++)
	  {
	    fractal_memory.level=level;
	    for(Group* pgroup : fractal_memory.all_groups[level])
	      force_shear_at_point(*pgroup,fractal);
	  }
	force_shear_at_particle(fractal_memory,fractal);
      }
    groups_level(fractal,fractal_memory.all_groups);
    fractal.timing(-1,44);
    if(fractal_memory.steps % fractal_memory.number_steps_out == 0)
      tree_dump(fractal_memory);
    fractal.timing(1,44);
    fractal.timing(-1,26);
    clean_up(fractal_memory,misc);
    // clean_up(fractal_memory,misc,fractal_ghost);
    fractal.timing(1,26);
    Full_Stop(fractal_memory,38);
    fractal.timing_lev(-1,0);
    fractal.timing(-1,28);
    gather_particles(fractal_memory,fractal);
    fractal.timing(1,28);
    fractal.timing_lev(1,0);
    if(fractal.get_halo_fixed()) 
      {
	fractal.timing(-1,25);
	halo_force_fixed(fractal);
	fractal.timing(1,25);
      }
    Full_Stop(fractal_memory,-1);
    FileFractal << "number of everything exiting Fractal "  << " " << Group::number_groups << " " << Point::number_points << "\n";
    FileFractal << " Total number of particles exiting Fractal " << Particle::number_particles << "\n";
    FileFractal << " Made It fractal_force " << fractal_memory.steps << " " << fractal_memory.p_mess->Clock()-fractal_memory.p_mess->WallTime << "\n";
    fractal_memory.p_file->FlushAll();
    if(fractal_memory.p_mess->FractalRank == 0)
      cerr << " Finished FractalForce " << fractal_memory.steps << "\n";
  }
  void Full_Stop(Fractal_Memory& mem,int number)
  {
    Full_Stop(mem,mem.p_mess->FractalWorld,number);
  }
  void Full_Stop(Fractal_Memory& mem,MPI_Comm& Comm,int number)
  {
    if(number >= 0)
      mem.p_fractal->timing(-1,number);
    mem.p_mess->Full_Stop_Do_Not_Argue(Comm);
    if(number >= 0)
      mem.p_fractal->timing(1,number);
  }
}
