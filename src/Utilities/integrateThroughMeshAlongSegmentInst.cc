//------------------------------------------------------------------------------
// Explicit instantiations.
//------------------------------------------------------------------------------
#include "integrateThroughMeshAlongSegment.cc"

namespace Spheral {

  template Dim<1>::Scalar integrateThroughMeshAlongSegment<Dim<1>, Dim<1>::Scalar>(const vector<vector<Dim<1>::Scalar> >& values, const Dim<1>::Vector& xmin, const Dim<1>::Vector& xmax, const vector<unsigned>& ncells, const Dim<1>::Vector& s0, const Dim<1>::Vector& s1);
  template Dim<2>::Scalar integrateThroughMeshAlongSegment<Dim<2>, Dim<2>::Scalar>(const vector<vector<Dim<2>::Scalar> >& values, const Dim<2>::Vector& xmin, const Dim<2>::Vector& xmax, const vector<unsigned>& ncells, const Dim<2>::Vector& s0, const Dim<2>::Vector& s1);
  template Dim<3>::Scalar integrateThroughMeshAlongSegment<Dim<3>, Dim<3>::Scalar>(const vector<vector<Dim<3>::Scalar> >& values, const Dim<3>::Vector& xmin, const Dim<3>::Vector& xmax, const vector<unsigned>& ncells, const Dim<3>::Vector& s0, const Dim<3>::Vector& s1);

}
