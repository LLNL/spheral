//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "BPLWraps/SuffixTraits.hh"
#include "BPLWraps/DataNameTraits.hh"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template<> const string DataNameTraits<Dim<1>::Scalar>::DataName = "Scalar";
  template<> const string DataNameTraits<Dim<1>::Vector>::DataName = "Vector";
  template<> const string DataNameTraits<Dim<1>::Tensor>::DataName = "Tensor";
  template<> const string DataNameTraits<Dim<1>::SymTensor>::DataName = "SymTensor";

  template<> const string DataNameTraits<Dim<2>::Vector>::DataName = "Vector";
  template<> const string DataNameTraits<Dim<2>::Tensor>::DataName = "Tensor";
  template<> const string DataNameTraits<Dim<2>::SymTensor>::DataName = "SymTensor";

  template<> const string DataNameTraits<Dim<3>::Vector>::DataName = "Vector";
  template<> const string DataNameTraits<Dim<3>::Tensor>::DataName = "Tensor";
  template<> const string DataNameTraits<Dim<3>::SymTensor>::DataName = "SymTensor";
}
