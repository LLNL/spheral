# -*-makefile-*-
# Spheral++ FractalGravity package makefile
#-------------------------------------------------------------------------------
PKGNAME = FractalGravity
PKGDIR = $(PKGNAME)/
LIBTARGET = lib$(PKGNAME).$(DYLIBEXT)
SRCTARGETS = \
	fftw++.cc \
	assign_density.cc \
	buffer_points.cc \
	check_for_edge_trouble.cc \
	check_high.cc \
	clean_up.cc\
	daughter_group_nina.cc \
	equivalence_class.cc \
	force_at_particle.cc \
	force_at_particle_sharp.cc \
	force_at_point.cc \
	force_at_point_initial.cc \
	fractal_gravity.cc \
	go_ahead.cc \
	heavies.cc\
	high_groups.cc \
	high_pairs.cc \
	high_points.cc \
	highest_level_group.cc \
	initial_forces.cc \
	initial_forces_sharp.cc \
	isolated_solver.cc \
	list_buffer.cc \
	neighbors_nina.cc \
	particle_lists.cc \
	periodic_solver.cc \
	poisson_solver.cc \
	potential_start.cc \
	power_spectrum.cc \
	sor_2.cc \
	sort_3.cc \
	timing.cc \
	tree_start.cc

LINKLIST = linklist
PYFFLELIST = pyfflelist
SPHERALTOP = ..
ARCH = $(SPHERALARCH)

EXTRAINCLUDES = -I.
EXTRALIBS = -lfftw3

#-------------------------------------------------------------------------------
include $(SPHERALTOP)/helpers/makefile_taumaster
include $(SPHERALTOP)/helpers/makefile_master
