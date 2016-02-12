# -*-makefile-*-
# Spheral++ FractalGravity package makefile
#-------------------------------------------------------------------------------
PKGNAME = FractalGravity
PKGDIR = $(PKGNAME)/
LIBTARGET = lib$(PKGNAME).$(DYLIBEXT)
F77TARGETS = \
	add_dens.f \
	assign_density.f \
	backwards.f \
	better_order.f \
	buffer_points.f \
	check_for_edge_trouble.f \
	check_high.f \
	compare.f \
	daughter_group_erika.f \
	dump_points.f \
	emergency_search.f \
	energies.f \
	equivalence_class.f \
	find_chains.f \
	find_go_ahead.f \
	find_empties.f \
	find_real_empties.f \
	find_smart_empties.f \
	find_high_pairs.f \
	find_high_points.f \
	find_neighbors.f \
	find_nina.f \
	find_possibles.f \
	find_uppers.f \
	find_uppers_erika.f \
	find_downers_nina.f \
	force_at_particle.f \
	force_at_point.f \
	fourn.f \
	fractal_gravity.f \
	high_groups.f \
	highest_level_group.f \
	isolated_solver.f \
	it_is_inside.f \
	just_checking.f \
	list_buffer.f \
	look_where.f \
	make_decisions_down.f \
	make_decisions_up.f \
	make_decisions.f \
	make_decisions_erika.f \
	make_sweep_list.f \
	negatives_groups.f \
	negatives_particles.f \
	negatives_points.f \
	offset.f \
	periodic_solver.f \
	poisson_solver.f \
	potential_start.f \
	remove_duplicates.f \
	rlft3.f \
	run_fractal.f \
	set_constant_real.f \
	slow_search.f \
	sor.f \
	sum_prod.f \
	sum_up.f \
	sweep_left.f \
	take_a_step.f \
	test_big.f \
	test_it_all.f \
	test_wrap.f \
	test_pos.f \
	timing.f \
	tree_start.f \
	wrapper.f

#	compare_sphere.f \
#	daughter_group.f \
#	daughter_group_nina.f \
#	daughter_group_exp.f \
#	daughter_group_hope.f \
#	find_neighbors_nina.f \
#	find_neighbors_1.f \
#	emergency_search_f.f \
#	find_uppers_erika_1.f \
#	find_downers_nina_1.f \
#	find_smart_uppers.f \
#	find_uppers_hope.f \
#	look_where_erika.f \
#	find_downers_exp.f \
#	make_decisions_0.f \
#	make_decisions_nina.f \
#	daughter_smart.f \
#	isolated_solver_dummy.f \
#	periodic_solver_dummy.f \
#	poisson_solver_recent.f \
#	fake.f \
#	sor_recent.f \
#	sor_2.f \

LINKLIST = linklist
PYFFLELIST = pyfflelist
SPHERALTOP = ..
ARCH = $(SPHERALARCH)

##EXTRAFLAGS = -DDEBUG
#-------------------------------------------------------------------------------
include $(SPHERALTOP)/helpers/makefile_taumaster
include $(SPHERALTOP)/helpers/makefile_master
