//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "weibullFlawDistribution.cc"

template Spheral::FieldSpace::Field<Spheral::Dim<1>, std::vector<double> > 
Spheral::PhysicsSpace::weibullFlawDistributionBenzAsphaug<Spheral::Dim<1> >(double,
                                                                            const double,
                                                                            const unsigned,
                                                                            const double,
                                                                            const double,
                                                                            const Spheral::NodeSpace::FluidNodeList<Spheral::Dim<1> >&,
                                                                            const int,
                                                                            const int);
template Spheral::FieldSpace::Field<Spheral::Dim<2>, std::vector<double> > 
Spheral::PhysicsSpace::weibullFlawDistributionBenzAsphaug<Spheral::Dim<2> >(double,
                                                                            const double,
                                                                            const unsigned,
                                                                            const double,
                                                                            const double,
                                                                            const Spheral::NodeSpace::FluidNodeList<Spheral::Dim<2> >&,
                                                                            const int,
                                                                            const int);
template Spheral::FieldSpace::Field<Spheral::Dim<3>, std::vector<double> > 
Spheral::PhysicsSpace::weibullFlawDistributionBenzAsphaug<Spheral::Dim<3> >(double,
                                                                            const double,
                                                                            const unsigned,
                                                                            const double,
                                                                            const double,
                                                                            const Spheral::NodeSpace::FluidNodeList<Spheral::Dim<3> >&,
                                                                            const int,
                                                                            const int);

template Spheral::FieldSpace::Field<Spheral::Dim<1>, std::vector<double> > 
Spheral::PhysicsSpace::weibullFlawDistributionOwen<Spheral::Dim<1> >(const unsigned,
                                                                     const double,
                                                                     const double,
                                                                     const Spheral::NodeSpace::FluidNodeList<Spheral::Dim<1> >&,
                                                                     const int,
                                                                     const double);
template Spheral::FieldSpace::Field<Spheral::Dim<2>, std::vector<double> > 
Spheral::PhysicsSpace::weibullFlawDistributionOwen<Spheral::Dim<2> >(const unsigned,
                                                                     const double,
                                                                     const double,
                                                                     const Spheral::NodeSpace::FluidNodeList<Spheral::Dim<2> >&,
                                                                     const int,
                                                                     const double);
template Spheral::FieldSpace::Field<Spheral::Dim<3>, std::vector<double> > 
Spheral::PhysicsSpace::weibullFlawDistributionOwen<Spheral::Dim<3> >(const unsigned,
                                                                     const double,
                                                                     const double,
                                                                     const Spheral::NodeSpace::FluidNodeList<Spheral::Dim<3> >&,
                                                                     const int,
                                                                     const double);
