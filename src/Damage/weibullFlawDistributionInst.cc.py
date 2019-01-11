text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "weibullFlawDistribution.cc"

template Spheral::Field<Spheral::Dim< %(ndim)s >, std::vector<double> > 
Spheral::weibullFlawDistributionBenzAsphaug<Spheral::Dim< %(ndim)s > >(double,
                                                                       const double,
                                                                       const unsigned,
                                                                       const double,
                                                                       const double,
                                                                       const Spheral::FluidNodeList<Spheral::Dim< %(ndim)s > >&,
                                                                       const int,
                                                                       const int);

template Spheral::Field<Spheral::Dim< %(ndim)s >, std::vector<double> > 
Spheral::weibullFlawDistributionOwen<Spheral::Dim< %(ndim)s > >(const unsigned,
                                                                const double,
                                                                const double,
                                                                const Spheral::FluidNodeList<Spheral::Dim< %(ndim)s > >&,
                                                                const int,
                                                                const double);
"""
