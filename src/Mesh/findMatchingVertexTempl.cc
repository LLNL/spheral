#include "findMatchingVertex.hh"

namespace Spheral {


//------------------------------------------------------------------------------
// Find the closest vertex in the list to the given position.
//------------------------------------------------------------------------------
template<typename Vector>
unsigned
findMatchingVertex(const Vector& target,
                   const std::vector<Vector>& verticesj) {
  const unsigned n = verticesj.size();
  unsigned i, result = n + 1;
  double chi2, chi2_min = 1e100;
  for (i = 0; i != n; ++i) {
    chi2 = (verticesj[i] - target).magnitude2();
    if (chi2 < chi2_min) {
      chi2_min = chi2;
      result = i;
    }
  }
  ENSURE(result < n);
  return result;
}

//------------------------------------------------------------------------------
// Find the closest vertex in a subset of the list to the given position.
//------------------------------------------------------------------------------
template<typename Vector>
unsigned
findMatchingVertex(const Vector& target,
                   const std::vector<Vector>& verticesj,
                   const std::vector<unsigned>& indicesj) {
  REQUIRE(indicesj.size() > 0);
  const unsigned n = indicesj.size();
  unsigned i, result = n + 1;
  double chi2, chi2_min = 1e100;
  for (i = 0; i != n; ++i) {
    CHECK(indicesj[i] < verticesj.size());
    chi2 = (verticesj[indicesj[i]] - target).magnitude2();
    if (chi2 < chi2_min) {
      chi2_min = chi2;
      result = i;
    }
  }
  ENSURE(result < n);
  return indicesj[result];
}

}
