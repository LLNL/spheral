//---------------------------------Spheral++----------------------------------//
// medianPosition
//
// Compute a definition of the median position for a collection of Vectors.
//
// Created by JMO, Thu Feb 18 11:25:55 PST 2010
//----------------------------------------------------------------------------//
#include <algorithm>
#include <vector>
#include <algorithm>

#include "medianPosition.hh"
#include "Geometry/Dimension.hh"

namespace Spheral {

using namespace std;

//------------------------------------------------------------------------------
// Compare two vectors by the given index coordinate.
//------------------------------------------------------------------------------
template<typename Vector>
struct CompareVectorsByIndex {
  CompareVectorsByIndex(const size_t index): mIndex(index) {}
  bool operator()(const Vector& lhs, const Vector& rhs) const {
    return lhs(mIndex) < rhs(mIndex);
  }
  size_t mIndex;
};

//------------------------------------------------------------------------------
// medianPosition
// Note this methd *will* shuffle the order of the input vector!
//------------------------------------------------------------------------------
template<typename Vector>
Vector
medianPosition(vector<Vector>& positions) {
  vector<Vector> p(positions);
  Vector result;
  const size_t n = p.size()/2;
  for (size_t idim = 0; idim != Vector::nDimensions; ++idim) {
    nth_element(p.begin(), p.begin() + n, p.end(), CompareVectorsByIndex<Vector>(idim));
    result(idim) = (*(p.begin() + n))(idim);
  }
  return result;
}

}
