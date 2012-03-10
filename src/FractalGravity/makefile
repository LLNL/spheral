# -*-makefile-*-
# Spheral++ FractalGravity package makefile
#-------------------------------------------------------------------------------
PKGNAME = FractalGravity
PKGDIR = $(PKGNAME)/
LIBTARGET = lib$(PKGNAME).$(DYLIBEXT)
SRCTARGETS = \
	FractalGravity.cc \
	adj_nina.cc \
	assign_density.cc \
	buffer_points.cc \
	candidate_points.cc \
	check_for_edge_trouble.cc \
	check_high.cc \
	clean_up.cc \
	cosmos.cc \
	daughter_group_nina.cc \
	equivalence_class.cc \
	fix_memory.cc \
	force_at_particle.cc \
	force_at_particle_other.cc \
	force_at_particle_sharp.cc \
	force_at_point.cc \
	force_at_point_initial.cc \
	force_shear_at_particle.cc \
	force_shear_at_point.cc \
	force_sum_particle.cc \
	force_test.cc \
	fractal_gravity.cc \
	fractal_gravity_wrapper.cc \
	go_ahead_points.cc \
	heavies.cc\
	high_groups.cc \
	high_pairs.cc \
	high_points.cc \
	initial_forces_sharp.cc \
	isolated_solver.cc \
	list_buffer.cc \
	max_predict.cc \
	neighbor_easy.cc \
	particle_lists.cc \
	particle_lists_fixed.cc \
	periodic_solver.cc \
	poisson_solver.cc \
	potential_start.cc \
	power_spectrum.cc \
	sort3_list.cc \
	split_particle.cc \
	test_group.cc \
	test_tree.cc \
	tree_start.cc \
	try_harder.cc \
	update_rv.cc \
	write_rv.cc \
	cosmos_power.cc \
	energy_leapfrog.cc \
	step_leapfrog.cc \
	halo_force_sph_harm.cc \
	halo_force_flat.cc \
	fractal_memory_parameters_cosmo.cc \
	make_particles_cosmo.cc \
	sor_2.cc \
	velocities.cc

LINKLIST = linklist
PYFFLELIST = pyfflelist
SPHERALTOP = ..
ARCH = $(SPHERALARCH)

EXTRAINCLUDES =
EXTRALIBS = -lfftw3 -lfftw++

#-------------------------------------------------------------------------------
include $(SPHERALTOP)/helpers/makefile_taumaster
include $(SPHERALTOP)/helpers/makefile_master
