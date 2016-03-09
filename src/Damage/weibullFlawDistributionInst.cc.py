text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "weibullFlawDistribution.cc"

template Spheral::FieldSpace::Field<Spheral::Dim< %(ndim)s >, std::vector<double> > 
Spheral::PhysicsSpace::weibullFlawDistributionBenzAsphaug<Spheral::Dim< %(ndim)s > >(double,
                                                                                     const double,
                                                                                     const unsigned,
                                                                                     const double,
                                                                                     const double,
                                                                                     const Spheral::NodeSpace::FluidNodeList<Spheral::Dim< %(ndim)s > >&,
                                                                                     const int,
                                                                                     const int);

template Spheral::FieldSpace::Field<Spheral::Dim< %(ndim)s >, std::vector<double> > 
Spheral::PhysicsSpace::weibullFlawDistributionOwen<Spheral::Dim< %(ndim)s > >(const unsigned,
                                                                              const double,
                                                                              const double,
                                                                              const Spheral::NodeSpace::FluidNodeList<Spheral::Dim< %(ndim)s > >&,
                                                                              const int,
                                                                              const double);
"""
