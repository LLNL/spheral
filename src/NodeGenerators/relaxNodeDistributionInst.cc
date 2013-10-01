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
                                      const WeightingFunctor<Dim<2> >& weighting,
                                      const int maxIterations,
                                      const double tolerance);
}
