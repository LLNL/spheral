//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "relaxNodeDistribution.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template relaxNodeDistribution(DataBaseSpace::DataBase<Dim<1> >& dataBase,
                                 const Dim<1>::FacetedVolume& boundary,
                                 const std::vector<BoundarySpace::Boundary<Dim<1> >*>& boundaries,
                                 const KernelSpace::TableKernel<Dim<1> >& W,
                                 const NodeSpace::SmoothingScaleBase<Dim<1> >& smoothingScaleMethod,
                                 const WeightingFunctor& weighting,
                                 const int maxIterations,
                                 const double tolerance);
}
