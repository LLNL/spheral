//---------------------------------Spheral++------------------------------------
// Compute the hull for a given points neighbor set.
//------------------------------------------------------------------------------
#include "computeNeighborHull.hh"
#include "Field/Field.hh"
#include "Field/FieldList.hh"

namespace Spheral {

using namespace std;

using FieldSpace::Field;
using FieldSpace::FieldList;

template<typename Dimension>
typename Dimension::FacetedVolume
computeNeighborHull(const std::vector<std::vector<int> >& fullConnectivity,
                    const typename Dimension::Scalar etaCutoff,
                    const typename Dimension::Vector& ri,
                    const typename Dimension::SymTensor& Hi,
                    const FieldSpace::FieldList<Dimension, typename Dimension::Vector>& position) {

  // Pre-conditions.
  const size_t numNodeLists = fullConnectivity.size();
  REQUIRE(etaCutoff > 0.0);
  REQUIRE(position.size() == numNodeLists);

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Dimension::FacetedVolume FacetedVolume;

  // Find the neighbor positions.
  const double etaCutoff2 = etaCutoff*etaCutoff;
  vector<Vector> neighbors(1, Vector::zero);
  for (size_t nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    for (size_t nodeListj = 0; nodeListj != numNodeLists; ++nodeListj) {
      const vector<int>& connectivity = fullConnectivity[nodeListj];
      for (vector<int>::const_iterator jItr = connectivity.begin();
           jItr != connectivity.end();
           ++jItr) {
        const int j = *jItr;
        const Vector& rj = position(nodeListj, j);
        const Vector rji = rj - ri;
        const Scalar etai2 = (Hi*rji).magnitude2();
        if (etai2 < etaCutoff2) neighbors.push_back(rji);
      }
    }
  }

  // Build the hull.
  return FacetedVolume(neighbors);
}

}

