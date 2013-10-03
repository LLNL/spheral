//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "relaxNodeDistribution.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template void relaxNodeDistribution(DataBaseSpace::DataBase<Dim<2> >& dataBase,
                                      const Dim<2>::FacetedVolume& boundary,
                                      const std::vector<BoundarySpace::Boundary<Dim<2> >*>& boundaries,
                                      const KernelSpace::TableKernel<Dim<2> >& W,
                                      const NodeSpace::SmoothingScaleBase<Dim<2> >& smoothingScaleMethod,
                                      const WeightingFunctor<Dim<2> >& weightingFunctor,
                                      const WeightingFunctor<Dim<2> >& massDensityFunctor,
                                      const double targetMass,
                                      const int maxIterations,
                                      const double tolerance);
  template void relaxNodeDistribution(DataBaseSpace::DataBase<Dim<3> >& dataBase,
                                      const Dim<3>::FacetedVolume& boundary,
                                      const std::vector<BoundarySpace::Boundary<Dim<3> >*>& boundaries,
                                      const KernelSpace::TableKernel<Dim<3> >& W,
                                      const NodeSpace::SmoothingScaleBase<Dim<3> >& smoothingScaleMethod,
                                      const WeightingFunctor<Dim<3> >& weightingFunctor,
                                      const WeightingFunctor<Dim<3> >& massDensityFunctor,
                                      const double targetMass,
                                      const int maxIterations,
                                      const double tolerance);
}
