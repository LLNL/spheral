//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "DEM/IncrementPairFieldList.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class IncrementPairFieldList<Dim<1>, std::vector<int>>;
  template class IncrementPairFieldList<Dim<1>, std::vector<Dim<1>::Scalar>>;
  template class IncrementPairFieldList<Dim<1>, std::vector<Dim<1>::Vector>>;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class IncrementPairFieldList<Dim<2>, std::vector<int>>;
  template class IncrementPairFieldList<Dim<2>, std::vector<Dim<2>::Scalar>>;
  template class IncrementPairFieldList<Dim<2>, std::vector<Dim<2>::Vector>>;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class IncrementPairFieldList<Dim<3>, std::vector<int>>;
  template class IncrementPairFieldList<Dim<3>, std::vector<Dim<3>::Scalar>>;
  template class IncrementPairFieldList<Dim<3>, std::vector<Dim<3>::Vector>>;
#endif
}