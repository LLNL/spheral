text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "Damage/weibullFlawDistribution.cc"

template Spheral::Field<Spheral::Dim< %(ndim)s >, std::vector<double> > 
Spheral::weibullFlawDistributionBenzAsphaug<Spheral::Dim< %(ndim)s > >(double,
                                                                       const double,
                                                                       const unsigned,
                                                                       const double,
                                                                       const double,
                                                                       const Spheral::FluidNodeList<Spheral::Dim< %(ndim)s > >&,
                                                                       const int,
                                                                       const int,
                                                                       const Spheral::Field<Spheral::Dim<%(ndim)s>, int>&);

template Spheral::Field<Spheral::Dim< %(ndim)s >, std::vector<double> > 
Spheral::weibullFlawDistributionOwen<Spheral::Dim< %(ndim)s > >(const unsigned,
                                                                const double,
                                                                const double,
                                                                const Spheral::FluidNodeList<Spheral::Dim< %(ndim)s > >&,
                                                                const int,
                                                                const double,
                                                                const Spheral::Field<Spheral::Dim<%(ndim)s>, int>&);
"""
