//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "NodeGenerators/relaxNodeDistribution.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template void relaxNodeDistribution(DataBase<Dim<1> >& dataBase,
                                      const Dim<1>::FacetedVolume& boundary,
                                      const std::vector<Boundary<Dim<1> >*>& boundaries,
                                      const TableKernel<Dim<1> >& W,
                                      const SmoothingScaleBase<Dim<1> >& smoothingScaleMethod,
                                      const WeightingFunctor<Dim<1> >& weightingFunctor,
                                      const WeightingFunctor<Dim<1> >& massDensityFunctor,
                                      const double targetMass,
                                      const int maxIterations,
                                      const double tolerance);
#endif

#if defined(SPHERAL_ENABLE_2D)
  template void relaxNodeDistribution(DataBase<Dim<2> >& dataBase,
                                      const Dim<2>::FacetedVolume& boundary,
                                      const std::vector<Boundary<Dim<2> >*>& boundaries,
                                      const TableKernel<Dim<2> >& W,
                                      const SmoothingScaleBase<Dim<2> >& smoothingScaleMethod,
                                      const WeightingFunctor<Dim<2> >& weightingFunctor,
                                      const WeightingFunctor<Dim<2> >& massDensityFunctor,
                                      const double targetMass,
                                      const int maxIterations,
                                      const double tolerance);
#endif

#if defined(SPHERAL_ENABLE_3D)
  template void relaxNodeDistribution(DataBase<Dim<3> >& dataBase,
                                      const Dim<3>::FacetedVolume& boundary,
                                      const std::vector<Boundary<Dim<3> >*>& boundaries,
                                      const TableKernel<Dim<3> >& W,
                                      const SmoothingScaleBase<Dim<3> >& smoothingScaleMethod,
                                      const WeightingFunctor<Dim<3> >& weightingFunctor,
                                      const WeightingFunctor<Dim<3> >& massDensityFunctor,
                                      const double targetMass,
                                      const int maxIterations,
                                      const double tolerance);
#endif
}