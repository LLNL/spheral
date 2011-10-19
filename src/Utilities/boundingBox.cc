//---------------------------------Spheral++----------------------------------//
// boundingBox
//
// Compute the minimum bounding box for a set of points.
//
// Created by JMO, Wed Oct 19 09:58:58 PDT 2011
//----------------------------------------------------------------------------//
#include <vector>
#include <algorithm>

#include "boundingBox.hh"
#include "Geometry/Dimension.hh"

namespace Spheral {

using namespace std;

//------------------------------------------------------------------------------
// Compute the minimum volume box containing all the points in a vector of
// positions.
// Note no global operations here!
//------------------------------------------------------------------------------
template<typename Vector>
void
boundingBox(const vector<Vector>& positions,
            Vector& xmin,
            Vector& xmax,
            const bool quantize) {

  xmin = DBL_MAX;
  xmax = -DBL_MAX;
  const unsigned n = positions.size();
  for (unsigned i = 0; i != n; ++i) {
    xmin = elementWiseMin(xmin, positions[i]);
    xmax = elementWiseMax(xmax, positions[i]);
  }

  // We make things integer values in an effort to make our result domain 
  // decomposition independent.
  if (quantize) {
    for (unsigned k = 0; k != Vector::nDimensions; ++k) {
      xmin(k) = double(int(xmin(k)) - (xmin(k) < 0.0 ? 1 : 0));
      xmax(k) = double(int(xmax(k)) + (xmax(k) < 0.0 ? 0 : 1));
    }
  }
}

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
template void boundingBox(const vector<Dim<1>::Vector>& positions,
                          Dim<1>::Vector& xmin,
                          Dim<1>::Vector& xmax,
                          const bool quantize);
template void boundingBox(const vector<Dim<2>::Vector>& positions,
                          Dim<2>::Vector& xmin,
                          Dim<2>::Vector& xmax,
                          const bool quantize);
template void boundingBox(const vector<Dim<3>::Vector>& positions,
                          Dim<3>::Vector& xmin,
                          Dim<3>::Vector& xmax,
                          const bool quantize);

}
