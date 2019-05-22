classes = fractal_memory.hh point_class.hh group_class.hh fractal_class.hh misc_class.hh particle_class.hh mess.hh file_class.hh kdtree.hh point_comp.hh fractal_interface_public.hh libs.hh headers.hh classes.hh
fopt1 = -pedantic -O3 
Wopt1 = -Wunused
gprof = 
debug = 
verbose = -v
debug = 	
verbose = 
STD = -std=c++11
dir   = execs/
MP = -qopenmp
#
system = quartz
fvers = 4
hvers = 11
#
compiler = mpiicpc
PREFIXF = /g/g90/jensv/cplus/fftw-3.3.4/install-intel
PREFIXH = /g/g90/jensv/cplus/HYPREINTEL110/hypre-2.11.0/src/hypre
#
ifeq ($(debug),d)
	fopt1 = -g -pedantic -O0
endif
#
ifeq ($(fvers),3)
	PREFIXF = /g/g90/jensv/cplus/fftw-3.3.3/install-intel
endif
#
ifeq ($(hvers),9)
	PREFIXH = /g/g90/jensv/cplus/hypre-2.9.0b/src/hypre/intel-global
endif
ifeq ($(hvers),10)
	PREFIXH = /g/g90/jensv/cplus/HYPREINTEL101/hypre-2.10.1/src/hypre
endif
#
ifeq ($(system),ibmgnu)
	MP = -fopenmp
	compiler =  mpig++-4.7.2
	PREFIXF = /g/g90/jensv/cplus/fftw-3.3.4/install-ibmgnu
	ifeq ($(fvers),3)
		PREFIXF = /g/g90/jensv/cplus/fftw-3.3.3/install-ibmgnu
	endif
	PREFIXH = /g/g90/jensv/cplus/HYPREIBM110/hypre-2.11.0/src/hypre
	ifeq ($(hvers),9)
		PREFIXH = /g/g90/jensv/cplus/hypreibmgnu/hypre-2.9.0b/src/hypre/noglobal
	endif
	ifeq ($(hvers),10)
		PREFIXH = /g/g90/jensv/cplus/HYPREIBM101/hypre-2.10.1/src/hypre
	endif
endif
ifeq ($(system),quartz)
	MP = -fopenmp
	compiler =  mpicc
	PREFIXF = /g/g90/jensv/cplus/fftw-3.3.7/intel-quartz
	PREFIXH = /g/g90/jensv/cplus/hypre-2.11.2/quartz/src/hypre
endif
Loptfftw = -L$(PREFIXF)/lib
head_dir_fftw = $(PREFIXF)/include
Lopthypre = -L$(PREFIXH)/lib
head_dir_hypre = $(PREFIXH)/include
#
compile = $(compiler) -I $(head_dir_fftw) -I $(head_dir_hypre) $(MP) $(STD) $(verbose) $(fopt1) $(Wopt1) $(gprof)
#
parameter_per_file = fractal_memory_parameters_cosmo
fractal_galaxy_run_o = fractal_galaxy.o
fractal_galaxy_spheral_o = fractal_galaxy_spheral.o
#
loptfftw = -lfftw3_mpi -lfftw3 -lm
lopthypre = -lHYPRE 
#lopthypre = $(Lopthypre) -lHYPRE 
#
wrapper_file = fractal_force_wrapper
cosmos_power_file = cosmos_power
energy_file = energy_leapfrog
halo_fixed_file = halo_force_flat
parameter_gal_file = fractal_memory_parameters_gal
parameter_isol_file = fractal_memory_parameters_gal
particle_gal_file = make_particles_gal
particle_cosmo_file = make_particles_cosmo
particle_cosmo_mpi_file = make_particles_cosmo_mpi
step_file = step_leapfrog
velocity_file = velocities
#
fractal_class.o :	fractal_class.cc $(classes)
	$(compile) -c fractal_class.cc
#
fractal_memory.o :	fractal_memory.cc $(classes)
	$(compile) -c fractal_memory.cc
#
file_class.o :	file_class.cc $(classes)
	$(compile) -c file_class.cc
#
group_class.o :	group_class.cc $(classes)
	$(compile) -c group_class.cc
#
mess.o :	mess.cc $(classes)
	$(compile) -c mess.cc
#
messA.o :	messA.cc $(classes)
	$(compile) -c messA.cc
#
misc_class.o :	misc_class.cc $(classes)
	$(compile) -c misc_class.cc
#
octree.o :	octree.cc $(classes)
	$(compile) -c octree.cc
#
kdtree.o :	kdtree.cc $(classes)
	$(compile) -c kdtree.cc
#
kdtreeA.o :	kdtreeA.cc $(classes)
	$(compile) -c kdtreeA.cc
#
particle_class.o :	particle_class.cc $(classes)
	$(compile) -c particle_class.cc
#
point_class.o :	point_class.cc $(classes)
	$(compile) -c point_class.cc
#
add_buffer_values.o :	add_buffer_values.cc $(classes)
	$(compile) -c add_buffer_values.cc
#
add_pseudo_particles.o :	add_pseudo_particles.cc $(classes)
	$(compile) -c add_pseudo_particles.cc
#
any_overlaps.o :	any_overlaps.cc $(classes)
	$(compile) -c any_overlaps.cc
#
adj_nina.o :	adj_nina.cc $(classes)
	$(compile) -c adj_nina.cc
#
am_I_conservative_enough_cosmo.o :	am_I_conservative_enough_cosmo.cc $(classes)
	$(compile) -c am_I_conservative_enough_cosmo.cc
#
am_I_conservative_enough_isol.o :	am_I_conservative_enough_isol.cc $(classes)
	$(compile) -c am_I_conservative_enough_isol.cc
#
assign_density.o :	assign_density.cc $(classes)
	$(compile) -c assign_density.cc
#
assign_densityB.o :	assign_densityB.cc $(classes)
	$(compile) -c assign_densityB.cc
#
balance_by_particles.o :	balance_by_particles.cc $(classes)
	$(compile) -c balance_by_particles.cc
#
balance_by_particles_cosmo.o :	balance_by_particles_cosmo.cc $(classes)
	$(compile) -c balance_by_particles_cosmo.cc
#
binary_balancing.o :	binary_balancing.cc $(classes)
	$(compile) -c binary_balancing.cc
#
box_stats.o :	box_stats.cc $(classes)
	$(compile) -c box_stats.cc
#
buffer_points.o :	buffer_points.cc $(classes)
	$(compile) -c buffer_points.cc
#
candidate_points.o :	candidate_points.cc $(classes)
	$(compile) -c candidate_points.cc
#
check_for_edge_trouble.o :	check_for_edge_trouble.cc $(classes)
	$(compile) -c check_for_edge_trouble.cc
#
check_high.o :	check_high.cc $(classes)
	$(compile) -c check_high.cc
#
clean_overlaps.o :	clean_overlaps.cc $(classes)
	$(compile) -c clean_overlaps.cc
#
clean_groups.o :	clean_groups.cc $(classes)
	$(compile) -c clean_groups.cc
#
clean_up.o :	clean_up.cc $(classes)
	$(compile) -c clean_up.cc
#
clean_container.o :	clean_container.cc $(classes)
	$(compile) -c clean_container.cc
#
cosmos.o :	cosmos.cc $(classes)
	$(compile) -c cosmos.cc
#
daughter_group_nina.o :	daughter_group_nina.cc $(classes)
	$(compile) -c daughter_group_nina.cc
#
dens_to_slices.o :	dens_to_slices.cc $(classes)
	$(compile) -c dens_to_slices.cc
#
dump_cosmo_boxes.o :	dump_cosmo_boxes.cc $(classes)
	$(compile) -c dump_cosmo_boxes.cc
#
dump_tree.o :	dump_tree.cc $(classes)
	$(compile) -c dump_tree.cc
#
energies.o :	energies.cc $(classes)
	$(compile) -c energies.cc
#
equivalence_class.o :	equivalence_class.cc $(classes)
	$(compile) -c equivalence_class.cc
#
fftw_where.o :	fftw_where.cc $(classes)
	$(compile) -c fftw_where.cc
#
find_global_level_max.o :	find_global_level_max.cc $(classes)
	$(compile) -c find_global_level_max.cc
#
findPointByPos.o :	findPointByPos.cc $(classes)
	$(compile) -c findPointByPos.cc
#
fix_memory.o :	fix_memory.cc $(classes)
	$(compile) -c fix_memory.cc
#
force_at_particle_other.o :	force_at_particle_other.cc $(classes)
	$(compile) -c force_at_particle_other.cc
#
force_at_particle.o :	force_at_particle.cc $(classes)
	$(compile) -c force_at_particle.cc
#
force_at_particle_sharp.o :	force_at_particle_sharp.cc $(classes)
	$(compile) -c force_at_particle_sharp.cc
#
force_at_point.o :	force_at_point.cc $(classes)
	$(compile) -c force_at_point.cc
#
force_at_point_initial.o :	force_at_point_initial.cc $(classes)
	$(compile) -c force_at_point_initial.cc
#
force_shear_at_particle.o :	force_shear_at_particle.cc $(classes)
	$(compile) -c force_shear_at_particle.cc
#
force_shear_at_point.o :	force_shear_at_point.cc $(classes)
	$(compile) -c force_shear_at_point.cc
#
force_sum_particle.o :	force_sum_particle.cc $(classes)
	$(compile) -c force_sum_particle.cc
#
force_test.o :	force_test.cc $(classes)
	$(compile) -c force_test.cc
#
fractal_cosmo.o :	fractal_cosmo.cc $(classes)
	$(compile) -c fractal_cosmo.cc
#
fractal_force.o :	fractal_force.cc  $(classes) 
	$(compile) -c fractal_force.cc
#
fractal_force_sync.o :	fractal_force_sync.cc  $(classes) 
	$(compile) -c fractal_force_sync.cc
#
fractal_galaxy_spheral.o :	fractal_galaxy_spheral.cc  $(classes) 
	$(compile) -c fractal_galaxy_spheral.cc
#
fractal_galaxy_intel.o :	fractal_galaxy_intel.cc  $(classes) 
	$(compile) -c fractal_galaxy_intel.cc
#
fractal_galaxy.o :	fractal_galaxy.cc  $(classes) 
	$(compile) -c fractal_galaxy.cc
#
FractalGravity.o :	FractalGravity.cc  $(classes) 
	$(compile) -c FractalGravity.cc
#
fractal_force_init.o :	fractal_force_init.cc  $(classes) 
	$(compile) -c fractal_force_init.cc
#
fractal_interface_public.o :	fractal_interface_public.cc $(classes)
	$(compile) -c fractal_interface_public.cc
#
fractal_nina.o :	fractal_nina.cc $(classes)
	$(compile) -c fractal_nina.cc
#
fractal_nina_cosmo.o :	fractal_nina_cosmo.cc $(classes)
	$(compile) -c fractal_nina_cosmo.cc
#
fractal_nina_galaxy.o :	fractal_nina_galaxy.cc $(classes)
	$(compile) -c fractal_nina_galaxy.cc
#
fractal_test_run.o :	fractal_test_run.cc  $(classes) 
	$(compile) -c fractal_test_run.cc
#
gather_particles.o :	gather_particles.cc $(classes)
	$(compile) -c gather_particles.cc
#
go_ahead_points.o :	go_ahead_points.cc $(classes)
	$(compile) -c go_ahead_points.cc
#
groups_level.o :	groups_level.cc $(classes)
	$(compile) -c groups_level.cc
#
high_pairs.o :	high_pairs.cc $(classes)
	$(compile) -c high_pairs.cc
#
high_points.o :	high_points.cc $(classes)
	$(compile) -c high_points.cc
#
heavies.o :	heavies.cc $(classes)
	$(compile) -c heavies.cc
#
high_groups.o :	high_groups.cc $(classes)
	$(compile) -c high_groups.cc
#
highest_level_group.o :	highest_level_group.cc $(classes)
	$(compile) -c highest_level_group.cc
#
highest_level_group_other.o :	highest_level_group_other.cc $(classes)
	$(compile) -c highest_level_group_other.cc
#
hypre_best_boxes.o :	hypre_best_boxes.cc $(classes)
	$(compile) -c hypre_best_boxes.cc
#
hypre_clever_boxes.o :	hypre_clever_boxes.cc $(classes)
	$(compile) -c hypre_clever_boxes.cc
#
hypre_dump.o :	hypre_dump.cc $(classes)
	$(compile) -c hypre_dump.cc
#
hypre_struct_load_balance.o :	hypre_struct_load_balance.cc $(classes)
	$(compile) -c hypre_struct_load_balance.cc
#
hypre_points_boxes.o :	hypre_points_boxes.cc $(classes)
	$(compile) -c hypre_points_boxes.cc
#
hypre_points_clean.o :	hypre_points_clean.cc $(classes)
	$(compile) -c hypre_points_clean.cc
#
hypre_points_zero.o :	hypre_points_zero.cc $(classes)
	$(compile) -c hypre_points_zero.cc
#
hypre_points_struct.o :	hypre_points_struct.cc $(classes)
	$(compile) -c hypre_points_struct.cc
#
hypre_solver.o :	hypre_solver.cc $(classes)
	$(compile) -c hypre_solver.cc
#
hypre_send_pots.o :	hypre_send_pots.cc $(classes)
	$(compile) -c hypre_send_pots.cc
#
hypre_solver_empty.o :	hypre_solver_empty.cc $(classes)
	$(compile) -c hypre_solver_empty.cc
#
hypre_solve_struct.o :	hypre_solve_struct.cc $(classes)
	$(compile) -c hypre_solve_struct.cc
#
hypre_solve_structA.o :	hypre_solve_structA.cc $(classes)
	$(compile) -c hypre_solve_structA.cc
#
hypre_test_boxes.o :	hypre_test_boxes.cc $(classes)
	$(compile) -c hypre_test_boxes.cc
#
hypre_world_create.o :	hypre_world_create.cc $(classes)
	$(compile) -c hypre_world_create.cc
#
hypre_world_destroy.o :	hypre_world_destroy.cc $(classes)
	$(compile) -c hypre_world_destroy.cc
#
info_to_slices.o :	info_to_slices.cc $(classes)
	$(compile) -c info_to_slices.cc
#
info_to_slices_to_pot_init.o :	info_to_slices_to_pot_init.cc $(classes)
	$(compile) -c info_to_slices_to_pot_init.cc
#
initial_forces_sharp.o :	initial_forces_sharp.cc $(classes)
	$(compile) -c initial_forces_sharp.cc
#
initial_forces_sharp_mpi.o :	initial_forces_sharp_mpi.cc $(classes)
	$(compile) -c initial_forces_sharp_mpi.cc
#
initial_forces_sharp_mpi_empty.o :	initial_forces_sharp_mpi_empty.cc $(classes)
	$(compile) -c initial_forces_sharp_mpi_empty.cc
#
isolated_solver.o :	isolated_solver.cc $(classes)
	$(compile) -c isolated_solver.cc
#
isolated_solver_mpi.o :	isolated_solver_mpi.cc $(classes)
	$(compile) -c isolated_solver_mpi.cc
#
it_is_inside.o :	it_is_inside.cc $(classes)
	$(compile) -c it_is_inside.cc
#
left_right.o :	left_right.cc $(classes)
	$(compile) -c left_right.cc 
#
list_buffer.o :	list_buffer.cc $(classes)
	$(compile) -c list_buffer.cc
#
make_me_a_box.o :	make_me_a_box.cc $(classes)
	$(compile) -c make_me_a_box.cc
#
make_me_a_galaxy.o :	make_me_a_galaxy.cc $(classes)
	$(compile) -c make_me_a_galaxy.cc
#
make_me_clumps.o :	make_me_clumps.cc $(classes)
	$(compile) -c make_me_clumps.cc
#
make_me_some_particles.o :	make_me_some_particles.cc $(classes)
	$(compile) -c make_me_some_particles.cc
#
make_sphere.o :	make_sphere.cc $(classes)
	$(compile) -c make_sphere.cc
#
match_edges.o :	match_edges.cc $(classes)
	$(compile) -c match_edges.cc
#
max_predict.o :	max_predict.cc $(classes)
	$(compile) -c max_predict.cc
#
mini_solve.o :	mini_solve.cc $(classes)
	$(compile) -c mini_solve.cc
#
move_small_boxes.o :	move_small_boxes.cc $(classes)
	$(compile) -c move_small_boxes.cc
#
neighbor_easy.o :	neighbor_easy.cc $(classes)
	$(compile) -c neighbor_easy.cc
#
neighbors_nina.o :	neighbors_nina.cc $(classes)
	$(compile) -c neighbors_nina.cc
#
node_groups_struct.o :	node_groups_struct.cc $(classes)
	$(compile) -c node_groups_struct.cc
#
on_edge.o :	on_edge.cc $(classes)
	$(compile) -c on_edge.cc
#
overlap.o :	overlap.cc $(classes)
	$(compile) -c overlap.cc
#
particle_lists.o :	particle_lists.cc $(classes)
	$(compile) -c particle_lists.cc
#
particle_lists_fixed.o :	particle_lists_fixed.cc $(classes)
	$(compile) -c particle_lists_fixed.cc
#
periodic_solver.o :	periodic_solver.cc $(classes)
	$(compile) -c periodic_solver.cc
#
periodic_solver_mpi.o :	periodic_solver_mpi.cc $(classes)
	$(compile) -c periodic_solver_mpi.cc
#
points_on_nodes.o :	points_on_nodes.cc $(classes)
	$(compile) -c points_on_nodes.cc
#
poisson_solver_struct.o :	poisson_solver_struct.cc $(classes)
	$(compile) -c poisson_solver_struct.cc
#
potential_start.o :	potential_start.cc $(classes)
	$(compile) -c potential_start.cc
#
power_spectrum.o :	power_spectrum.cc $(classes)
	$(compile) -c power_spectrum.cc
#
problem_points.o :	problem_points.cc $(classes)
	$(compile) -c problem_points.cc
#
really_clear.o :	really_clear.cc $(classes)
	$(compile) -c really_clear.cc
#
remove_dupe_points.o :	remove_dupe_points.cc $(classes)
	$(compile) -c remove_dupe_points.cc
#
remove_pseudo_particles.o :	remove_pseudo_particles.cc $(classes)
	$(compile) -c remove_pseudo_particles.cc
#
right_diff.o :	right_diff.cc $(classes)
	$(compile) -c right_diff.cc
#
scatter_particles.o :	scatter_particles.cc $(classes)
	$(compile) -c scatter_particles.cc
#
scatter_particles_faster.o :	scatter_particles_faster.cc $(classes)
	$(compile) -c scatter_particles_faster.cc
#
scatter_particles_find.o :	scatter_particles_find.cc $(classes)
	$(compile) -c scatter_particles_find.cc
#
scatter_particles_simple.o :	scatter_particles_simple.cc $(classes)
	$(compile) -c scatter_particles_simple.cc
#
shortest_vector.o :	shortest_vector.cc $(classes)
	$(compile) -c shortest_vector.cc
#
shrink_cube.o :	shrink_cube.cc $(classes)
	$(compile) -c shrink_cube.cc
#
small_exceptions.o :	small_exceptions.cc $(classes)
	$(compile) -c small_exceptions.cc
#
sor_solver.o :	sor_solver.cc $(classes)
	$(compile) -c sor_solver.cc
#
sor_2.o :	sor_2.cc $(classes)
	$(compile) -c sor_2.cc
#
sort_3.o :	sort_3.cc $(classes)
	$(compile) -c sort_3.cc
#
slices_to_potf.o :	slices_to_potf.cc $(classes)
	$(compile) -c slices_to_potf.cc
#
slices_to_pot_init.o :	slices_to_pot_init.cc $(classes)
	$(compile) -c slices_to_pot_init.cc
#
split_nodes.o :	split_nodes.cc $(classes)
	$(compile) -c split_nodes.cc
#
split_particle.o :	split_particle.cc $(classes)
	$(compile) -c split_particle.cc
#
sort3_list.o :	sort3_list.cc $(classes)
	$(compile) -c sort3_list.cc
#
start_writing.o :	start_writing.cc $(classes)
	$(compile) -c start_writing.cc
#
sum_pot_forces.o :	sum_pot_forces.cc $(classes)
	$(compile) -c sum_pot_forces.cc
#
super_groups.o :	super_groups.cc $(classes)
	$(compile) -c super_groups.cc
#
take_a_cosmic_leap.o :	take_a_cosmic_leap.cc $(classes)
	$(compile) -c take_a_cosmic_leap.cc
#
take_a_leap_isol.o :	take_a_leap_isol.cc $(classes)
	$(compile) -c take_a_leap_isol.cc
#
testdir.o :	testdir.cc 
	$(compile) -c testdir.cc
#
testfftw.o :	testfftw.cc 
	$(compile) -c testfftw.cc
#
test3dfftw.o :	test3dfftw.cc 
	$(compile) -c test3dfftw.cc
#
testfftwS.o :	testfftwS.cc 
	$(compile) -c testfftwS.cc
#
test_gal.o :	test_gal.cc $(classes)
	$(compile) -c test_gal.cc
#
test_good_point.o :	test_good_point.cc $(classes)
	$(compile) -c test_good_point.cc
#
test_group.o :	test_group.cc $(classes)
	$(compile) -c test_group.cc
#
test_points.o :	test_points.cc $(classes)
	$(compile) -c test_points.cc
#
test_tree.o :	test_tree.cc $(classes)
	$(compile) -c test_tree.cc
#
total_mass.o :	total_mass.cc $(classes)
	$(compile) -c total_mass.cc
#
tree_dump.o :	tree_dump.cc $(classes)
	$(compile) -c tree_dump.cc
#
tree_start.o :	tree_start.cc $(classes)
	$(compile) -c tree_start.cc
#
tree_start_mini.o :	tree_start_mini.cc $(classes)
	$(compile) -c tree_start_mini.cc
#
try_harder.o :	try_harder.cc $(classes)
	$(compile) -c try_harder.cc
#
update_rv.o :	update_rv.cc $(classes)
	$(compile) -c update_rv.cc
#
use_freenodes.o :	use_freenodes.cc $(classes)
	$(compile) -c use_freenodes.cc
#
which_element.o :	which_element.cc $(classes)
	$(compile) -c which_element.cc
#
wrap.o :	wrap.cc $(classes)
	$(compile) -c wrap.cc
#
write_rv.o :	write_rv.cc $(classes)
	$(compile) -c write_rv.cc
#
$(cosmos_power_file).o :	$(cosmos_power_file).cc $(classes)
	$(compile) -c $(cosmos_power_file).cc
#
$(energy_file).o :	$(energy_file).cc $(classes) 
	$(compile) -c $(energy_file).cc
#
$(wrapper_file).o :	$(wrapper_file).cc $(classes) 
	$(compile) -c $(wrapper_file).cc
#
$(halo_fixed_file).o :	$(halo_fixed_file).cc $(classes)
	$(compile) -c $(halo_fixed_file).cc
#
$(parameter_per_file).o :	$(parameter_per_file).cc $(classes)
	$(compile) -c $(parameter_per_file).cc
#
$(parameter_gal_file).o :	$(parameter_gal_file).cc $(classes)
	$(compile) -c $(parameter_gal_file).cc
#
$(particle_cosmo_file).o :	$(particle_cosmo_file).cc $(classes)
	$(compile) -c $(particle_cosmo_file).cc
#
$(particle_cosmo_mpi_file).o :	$(particle_cosmo_mpi_file).cc $(classes)
	$(compile) -c $(particle_cosmo_mpi_file).cc
#
$(particle_gal_file).o :	$(particle_gal_file).cc $(classes)
	$(compile) -c $(particle_gal_file).cc
#
$(step_file).o :	$(step_file).cc $(classes)
	$(compile) -c $(step_file).cc
#
$(velocity_file).o :	$(velocity_file).cc $(classes)
	$(compile) -c $(velocity_file).cc
#
files_class_o= \
file_class.o \
fractal_class.o \
fractal_memory.o \
group_class.o \
mess.o \
misc_class.o \
kdtree.o \
particle_class.o \
point_class.o 
#
files_fractal_nina_o= \
add_pseudo_particles.o \
adj_nina.o \
am_I_conservative_enough_isol.o \
assign_density.o \
balance_by_particles.o \
binary_balancing.o \
buffer_points.o \
candidate_points.o \
check_for_edge_trouble.o \
check_high.o \
clean_up.o \
clean_container.o \
cosmos.o \
daughter_group_nina.o \
dump_cosmo_boxes.o \
dump_tree.o \
equivalence_class.o \
find_global_level_max.o \
fix_memory.o \
force_at_particle.o \
force_at_particle_sharp.o \
force_at_point.o \
force_at_point_initial.o \
force_shear_at_particle.o \
force_shear_at_point.o \
force_sum_particle.o \
force_test.o \
fractal_force.o \
fractal_force_init.o \
fractal_interface_public.o \
go_ahead_points.o \
groups_level.o \
high_groups.o \
high_pairs.o \
high_points.o \
left_right.o \
list_buffer.o \
make_me_a_galaxy.o \
match_edges.o \
max_predict.o \
mini_solve.o \
move_small_boxes.o \
neighbor_easy.o \
on_edge.o \
overlap.o \
potential_start.o \
power_spectrum.o \
remove_pseudo_particles.o \
shrink_cube.o \
small_exceptions.o \
sort3_list.o \
split_particle.o \
start_writing.o \
sum_pot_forces.o \
super_groups.o \
take_a_leap_isol.o \
test_gal.o \
test_group.o \
test_points.o \
test_tree.o \
total_mass.o \
tree_dump.o \
tree_start.o \
tree_start_mini.o \
try_harder.o \
update_rv.o \
use_freenodes.o \
write_rv.o
#
files_hypre_nina_o= \
add_buffer_values.o \
any_overlaps.o \
box_stats.o \
clean_overlaps.o \
hypre_best_boxes.o \
hypre_points_boxes.o \
hypre_points_clean.o \
hypre_points_struct.o \
hypre_points_zero.o \
hypre_solve_struct.o \
hypre_struct_load_balance.o \
hypre_test_boxes.o \
hypre_world_create.o \
node_groups_struct.o \
points_on_nodes.o \
poisson_solver_struct.o \
#
files_hypre_empty_o= \
hypre_solver_empty.o
#
files_fftw_mpi_nina_o= \
initial_forces_sharp_mpi.o \
isolated_solver_mpi.o \
periodic_solver_mpi.o 
#
files_fftw_nina_o= \
initial_forces_sharp.o \
isolated_solver.o \
periodic_solver.o 
#
files_mpi_comm_nina_o= \
dens_to_slices.o \
gather_particles.o \
scatter_particles_find.o \
slices_to_potf.o 
#

fractal_cosmo : 	 $(files_hypre_nina_o) $(files_class_o) $(files_fractal_nina_o) $(files_fftw_mpi_nina_o) $(files_mpi_comm_nina_o) fractal_cosmo.o $(cosmos_power_file).o $(energy_file).o $(halo_fixed_file).o $(parameter_per_file).o $(particle_cosmo_mpi_file).o $(step_file).o $(velocity_file).o $(wrapper_file).o
	$(compile) $(Loptfftw) $(Lopthypre) -o $(dir)fractal_nina_cosmo.exe $(files_class_o) $(files_fractal_nina_o) $(files_hypre_nina_o) $(files_fftw_mpi_nina_o) $(files_mpi_comm_nina_o) fractal_cosmo.o $(cosmos_power_file).o $(energy_file).o $(halo_fixed_file).o $(parameter_per_file).o $(particle_cosmo_mpi_file).o $(step_file).o $(velocity_file).o $(wrapper_file).o  $(lopthypre) $(loptfftw)

#
fractal_galaxy :	$(files_hypre_nina_o) $(files_class_o) $(files_fractal_nina_o) $(files_fftw_mpi_nina_o) $(files_mpi_comm_nina_o) $(fractal_galaxy_run_o) $(cosmos_power_file).o $(halo_fixed_file).o
	$(compile) $(Loptfftw) $(Lopthypre) -o $(dir)fractal_nina_galaxy.exe $(files_hypre_nina_o) $(files_class_o) $(files_fractal_nina_o) $(files_fftw_mpi_nina_o) $(files_mpi_comm_nina_o) $(fractal_galaxy_run_o) $(cosmos_power_file).o $(halo_fixed_file).o $(lopthypre) $(loptfftw)

fractal_galaxy_spheral :	$(files_hypre_nina_o) $(files_class_o) $(files_fractal_nina_o) $(files_fftw_mpi_nina_o) $(files_mpi_comm_nina_o) $(fractal_galaxy_spheral_o) $(cosmos_power_file).o $(halo_fixed_file).o
	$(compile) $(Loptfftw) $(Lopthypre) -o $(dir)fractal_nina_galaxy_spheral.exe $(files_hypre_nina_o) $(files_class_o) $(files_fractal_nina_o) $(files_fftw_mpi_nina_o) $(files_mpi_comm_nina_o) $(fractal_galaxy_spheral_o) $(cosmos_power_file).o $(halo_fixed_file).o $(lopthypre) $(loptfftw)
