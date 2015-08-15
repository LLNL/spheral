text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "relaxNodeDistribution.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template void relaxNodeDistribution(DataBaseSpace::DataBase<Dim< %(ndim)s > >& dataBase,
                                      const Dim< %(ndim)s >::FacetedVolume& boundary,
                                      const std::vector<BoundarySpace::Boundary<Dim< %(ndim)s > >*>& boundaries,
                                      const KernelSpace::TableKernel<Dim< %(ndim)s > >& W,
                                      const NodeSpace::SmoothingScaleBase<Dim< %(ndim)s > >& smoothingScaleMethod,
                                      const WeightingFunctor<Dim< %(ndim)s > >& weightingFunctor,
                                      const WeightingFunctor<Dim< %(ndim)s > >& massDensityFunctor,
                                      const double targetMass,
                                      const int maxIterations,
                                      const double tolerance);
}
"""
