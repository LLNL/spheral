classes = fractal_memory.hh classes.hh point_class.hh group_class.hh fractal_class.hh misc_class.hh particle_class.hh mess.hh file_class.hh libs.hh headers.hh fractal_interface_public.hh mini_Point_class.hh
fopt1 = -pedantic -O3 
Wopt1 = -Wunused
fftopt1 = 
gprof = 
debug = -g
verbose = -v
debug= 
verbose = 
STD = -std=c++0x
dir   = execs/
compiler = mpiicpc
#
system = intel
manynodes = short
#
HypreLib = /g/g90/jensv/cplus/hypre-2.9.0b/src/hypre/intel-global/lib
HypreDir = /g/g90/jensv/cplus/hypre-2.9.0b/src/hypre/intel-global/include
#
ifeq ($(system),intel)
	compiler = mpiicpc
	Loptfftw = -L/g/g90/jensv/cplus/fftw-3.3.3/install-intel/lib
	head_dir_fftw = /g/g90/jensv/cplus/fftw-3.3.3/install-intel/include
	ifeq ($(manynodes),short)
		Lopthypre = -L/g/g90/jensv/cplus/hypre-2.9.0b/src/hypre/intel-global/lib
		head_dir_hypre = /g/g90/jensv/cplus/hypre-2.9.0b/src/hypre/intel-global/include
	else
		Lopthypre = -L/g/g90/jensv/cplus/hypre-2.9.0b/src/hypre/intel-noglobal/lib
		head_dir_hypre = /g/g90/jensv/cplus/hypre-2.9.0b/src/hypre/intel-noglobal/include
	endif
endif
#
ifeq ($(system),ibmgnu)
	compiler = mpicxx
	Loptfftw = -L/g/g90/jensv/cplus/fftw-3.3.3/install-ibmgnu/lib
	head_dir_fftw = /g/g90/jensv/cplus/fftw-3.3.3/install-ibmgnu/include
	ifeq ($(manynodes),short)
		Lopthypre = -L/g/g90/jensv/cplus/hypreibmgnu/hypre-2.9.0b/src/hypre/global/lib
		head_dir_hypre = /g/g90/jensv/cplus/hypreibmgnu/hypre-2.9.0b/src/hypre/global/include
	else
		Lopthypre = -L/g/g90/jensv/cplus/hypreibmgnu/hypre-2.9.0b/src/hypre/noglobal/lib
		head_dir_hypre = /g/g90/jensv/cplus/hypreibmgnu/hypre-2.9.0b/src/hypre/noglobal/include
	endif
endif
#
compile = $(compiler) -I $(head_dir_fftw) -I $(head_dir_hypre) $(STD) $(verbose) $(debug) $(fopt1) $(Wopt1) $(gprof)
#
parameter_per_file = fractal_memory_parameters_cosmo
fractal_galaxy_run_o=fractal_galaxy.o
#
loptfftw = -lfftw3_mpi -lfftw3 -lm
lopthypre = $(Lopthypre) -lHYPRE 
#
#compile = $(compiler) -I $(head_dir_fftw) -I $(head_dir_hypre) $(STD) $(verbose) $(debug)
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
add_pseudo_particles.o :	add_pseudo_particles.cc $(classes)
	$(compile) -c add_pseudo_particles.cc
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
balance_by_particles.o :	balance_by_particles.cc $(classes)
	$(compile) -c balance_by_particles.cc
#
balance_by_particles_iter.o :	balance_by_particles_iter.cc $(classes)
	$(compile) -c balance_by_particles_iter.cc
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
clean_up.o :	clean_up.cc $(classes)
	$(compile) -c clean_up.cc
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
dens_to_slicesA.o :	dens_to_slicesA.cc $(classes)
	$(compile) -c dens_to_slicesA.cc
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
fractal_force.o :	fractal_force.cc  $(classes) 
	$(compile) -c fractal_force.cc
#
fractal_force_sync.o :	fractal_force_sync.cc  $(classes) 
	$(compile) -c fractal_force_sync.cc
#
fractal_galaxy_ibm.o :	fractal_galaxy_ibm.cc  $(classes) 
	$(compile) -c fractal_galaxy_ibm.cc
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
fractal_nbody.o :	fractal_nbody.cc $(classes)
	$(compile) -c fractal_nbody.cc
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
hypre_dump.o :	hypre_dump.cc $(classes)
	$(compile) -c hypre_dump.cc
#
hypre_ij_numbering.o :	hypre_ij_numbering.cc $(classes)
	$(compile) -c hypre_ij_numbering.cc
#
hypre_ij_numbering_find.o :	hypre_ij_numbering_find.cc $(classes)
	$(compile) -c hypre_ij_numbering_find.cc
#
hypre_ij_numbering_findA.o :	hypre_ij_numbering_findA.cc $(classes)
	$(compile) -c hypre_ij_numbering_findA.cc
#
hypre_ij_numbering_selfie.o :	hypre_ij_numbering_selfie.cc $(classes)
	$(compile) -c hypre_ij_numbering_selfie.cc
#
hypre_ij_solver.o :	hypre_ij_solver.cc $(classes)
	$(compile) -c hypre_ij_solver.cc
#
hypre_ij_solverA.o :	hypre_ij_solverA.cc $(classes)
	$(compile) -c hypre_ij_solverA.cc
#
hypre_ij_solver_selfie.o :	hypre_ij_solver_selfie.cc $(classes)
	$(compile) -c hypre_ij_solver_selfie.cc
#
hypre_ij_solver_pcg.o :	hypre_ij_solver_pcg.cc $(classes)
	$(compile) -c hypre_ij_solver_pcg.cc
#
hypre_load_balance.o :	hypre_load_balance.cc $(classes)
	$(compile) -c hypre_load_balance.cc
#
hypre_load_balanceA.o :	hypre_load_balanceA.cc $(classes)
	$(compile) -c hypre_load_balanceA.cc
#
hypre_load_balanceB.o :	hypre_load_balanceB.cc $(classes)
	$(compile) -c hypre_load_balanceB.cc
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
hypre_struct_solver.o :	hypre_struct_solver.cc $(classes)
	$(compile) -c hypre_struct_solver.cc
#
info_to_slices.o :	info_to_slices.cc $(classes)
	$(compile) -c info_to_slices.cc
#
info_to_slicesA.o :	info_to_slicesA.cc $(classes)
	$(compile) -c info_to_slicesA.cc
#
initial_forces_sharp.o :	initial_forces_sharp.cc $(classes)
	$(compile) -c initial_forces_sharp.cc
#
initial_forces_sharp_mpi.o :	initial_forces_sharp_mpi.cc $(classes)
	$(compile) -c initial_forces_sharp_mpi.cc
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
make_me_a_galaxy.o :	make_me_a_galaxy.cc $(classes)
	$(compile) -c make_me_a_galaxy.cc
#
make_me_some_particles.o :	make_me_some_particles.cc $(classes)
	$(compile) -c make_me_some_particles.cc
#
make_sphere.o :	make_sphere.cc $(classes)
	$(compile) -c make_sphere.cc
#
max_predict.o :	max_predict.cc $(classes)
	$(compile) -c max_predict.cc
#
neighbor_easy.o :	neighbor_easy.cc $(classes)
	$(compile) -c neighbor_easy.cc
#
neighbors_nina.o :	neighbors_nina.cc $(classes)
	$(compile) -c neighbors_nina.cc
#
new_balance.o :	new_balance.cc $(classes)
	$(compile) -c new_balance.cc
#
new_balanceA.o :	new_balanceA.cc $(classes)
	$(compile) -c new_balanceA.cc
#
node_split.o :	node_split.cc $(classes)
	$(compile) -c node_split.cc
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
potential_start.o :	potential_start.cc $(classes)
	$(compile) -c potential_start.cc
#
poisson_solver.o :	poisson_solver.cc $(classes)
	$(compile) -c poisson_solver.cc
#
poisson_solverA.o :	poisson_solverA.cc $(classes)
	$(compile) -c poisson_solverA.cc
#
power_spectrum.o :	power_spectrum.cc $(classes)
	$(compile) -c power_spectrum.cc
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
shortest_vector.o :	shortest_vector.cc $(classes)
	$(compile) -c shortest_vector.cc
#
sor_solver.o :	sor_solver.cc $(classes)
	$(compile) -c sor_solver.cc
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
sor_2.o :	sor_2.cc $(classes)
	$(compile) -c sor_2.cc
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
test_group.o :	test_group.cc $(classes)
	$(compile) -c test_group.cc
#
test_tree.o :	test_tree.cc $(classes)
	$(compile) -c test_tree.cc
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
#
files_fractal_nina_o= \
add_pseudo_particles.o \
adj_nina.o \
am_I_conservative_enough_isol.o \
assign_density.o \
buffer_points.o \
candidate_points.o \
check_for_edge_trouble.o \
check_high.o \
clean_up.o \
cosmos.o \
daughter_group_nina.o \
dump_tree.o \
equivalence_class.o \
find_global_level_max.o \
findPointByPos.o \
fix_memory.o \
force_at_particle.o \
force_at_particle_other.o \
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
heavies.o \
high_groups.o \
high_pairs.o \
high_points.o \
left_right.o \
list_buffer.o \
make_me_a_galaxy.o \
max_predict.o \
neighbor_easy.o \
new_balanceA.o \
on_edge.o \
overlap.o \
particle_lists.o \
particle_lists_fixed.o \
poisson_solver.o \
potential_start.o \
power_spectrum.o \
remove_pseudo_particles.o \
sor_solver.o \
sor_2.o \
sort3_list.o \
split_nodes.o \
split_particle.o \
start_writing.o \
sum_pot_forces.o \
take_a_cosmic_leap.o \
take_a_leap_isol.o \
test_gal.o \
test_group.o \
test_tree.o \
tree_start.o \
tree_start_mini.o \
try_harder.o \
update_rv.o \
write_rv.o
#
files_hypre_nina_o= \
hypre_dump.o \
hypre_ij_numbering_find.o \
hypre_ij_numbering_selfie.o \
hypre_ij_solver.o \
hypre_ij_solver_selfie.o \
hypre_load_balanceB.o \
hypre_send_pots.o
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
info_to_slices.o \
scatter_particles_find.o \
slices_to_potf.o \
slices_to_pot_init.o 
#
fractal_cosmo : 	$(files_fractal_nina_o) $(files_hypre_nina_o) $(files_fftw_mpi_nina_o) $(files_mpi_comm_nina_o) fractal_nbody.o $(cosmos_power_file).o $(energy_file).o $(halo_fixed_file).o $(parameter_per_file).o $(particle_cosmo_mpi_file).o $(step_file).o $(velocity_file).o $(wrapper_file).o
	$(compile) $(Loptfftw) $(Lopthypre) -o $(dir)fractal_nina_cosmo.exe $(files_fractal_nina_o) $(files_hypre_nina_o) $(files_fftw_mpi_nina_o) $(files_mpi_comm_nina_o) fractal_nbody.o $(cosmos_power_file).o $(energy_file).o $(halo_fixed_file).o $(parameter_per_file).o $(particle_cosmo_mpi_file).o $(step_file).o $(velocity_file).o $(wrapper_file).o  $(lopthypre) $(loptfftw)
#
#fractal_per_single : 	$(files_fractal_nina_o) $(files_hypre_empty_o) $(files_fftw_nina_o) $(files_mpi_comm_nina_o) fractal_nbody.o $(cosmos_power_file).o $(energy_file).o $(halo_fixed_file).o $(parameter_per_file).o $(particle_cosmo_file).o $(step_file).o $(velocity_file).o $(wrapper_file).o
#	$(compile) $(Loptfftw) -o $(dir)fractal_nina_per_single.exe $(files_fractal_nina_o) $(files_hypre_empty_o) $(files_fftw_nina_o) $(files_mpi_comm_nina_o) fractal_nbody.o $(cosmos_power_file).o $(energy_file).o $(halo_fixed_file).o $(parameter_per_file).o $(particle_cosmo_file).o $(step_file).o $(velocity_file).o $(wrapper_file).o  $(loptfftw)
#
#fractal_gal_single : 	$(files_fractal_nina_o) $(files_hypre_empty_o) $(files_fftw_nina_o) $(files_mpi_comm_nina_o) fractal_nbody.o $(cosmos_power_file).o $(energy_file).o $(halo_fixed_file).o $(parameter_gal_file).o $(particle_gal_file).o $(step_file).o $(velocity_file).o $(wrapper_file).o
#	$(compile) $(Loptfftw) -o $(dir)fractal_nina_gal_single.exe $(files_fractal_nina_o) $(files_hypre_empty_o) $(files_fftw_nina_o) $(files_mpi_comm_nina_o) fractal_nbody.o $(cosmos_power_file).o $(energy_file).o $(halo_fixed_file).o $(parameter_gal_file).o $(particle_gal_file).o $(step_file).o $(velocity_file).o $(wrapper_file).o $(loptfftw)
#
#fractal_gal_mpi : 	$(files_fractal_nina_o) $(files_hypre_nina_o) $(files_fftw_mpi_nina_o) $(files_mpi_comm_nina_o) fractal_nbody.o $(cosmos_power_file).o $(energy_file).o $(halo_fixed_file).o $(parameter_gal_file).o $(particle_gal_file).o $(step_file).o $(velocity_file).o $(wrapper_file).o
#	$(compile) $(Loptfftw) $(Lopthypre) -o $(dir)fractal_nina_gal_mpi.exe $(files_fractal_nina_o) $(files_hypre_nina_o) $(files_fftw_mpi_nina_o) $(files_mpi_comm_nina_o) fractal_nbody.o $(cosmos_power_file).o $(energy_file).o $(halo_fixed_file).o $(parameter_gal_file).o $(particle_gal_file).o $(step_file).o $(velocity_file).o $(wrapper_file).o  $(lopthypre) $(loptfftw)
#
fractal_galaxy :	$(files_fractal_nina_o) $(files_hypre_nina_o) $(files_fftw_mpi_nina_o) $(files_mpi_comm_nina_o) $(fractal_galaxy_run_o) $(cosmos_power_file).o $(halo_fixed_file).o
	$(compile) $(Loptfftw) $(Lopthypre) -o $(dir)fractal_nina_galaxy.exe $(files_fractal_nina_o) $(files_hypre_nina_o) $(files_fftw_mpi_nina_o) $(files_mpi_comm_nina_o) $(fractal_galaxy_run_o) $(cosmos_power_file).o $(halo_fixed_file).o $(lopthypre) $(loptfftw)
#
#fractal_test_run_mpi :	$(files_fractal_nina_o) $(files_hypre_nina_o) $(files_fftw_mpi_nina_o) $(files_mpi_comm_nina_o) fractal_test_run.o $(cosmos_power_file).o $(halo_fixed_file).o
#	$(compile) $(Loptfftw) $(Lopthypre) -o $(dir)fractal_nina_test_mpi.exe $(files_fractal_nina_o) $(files_hypre_nina_o) $(files_fftw_mpi_nina_o) $(files_mpi_comm_nina_o) fractal_test_run.o $(cosmos_power_file).o $(halo_fixed_file).o $(lopthypre) $(loptfftw)
#
#
#fftopt1 = -DNDEBUG
#ATM Loptfftw = -L/g/g90/jensv/cplus/fftw-3.3.2/mpi/.libs
#Loptfftw = -L/usr/local/tools/fftw3-3.2/lib
#Lopthypre = -L/g/g90/jensv/cplus/hypre-2.9.0b/src/lib/
#fftwmpifile = /g/g90/jensv/cplus/fftw-3.3.2/mpi/.libs/libfftw3_mpi.a
#hyprempifile = /g/g90/jensv/cplus/hypre-2.8.0b/src/lib/libHYPRE.a
#ATM head_dir_fftw = /g/g90/jensv/cplus/fftw-3.3.2/mpi/
#head_dir_fftw = /usr/local/tools/fftw3-3.2/include
#head_dir_hypre = /g/g90/jensv/cplus/hypre-2.9.0b/src/hypre/include
#ATM compiler = mpiicpc -mpicpp
