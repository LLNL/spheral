//------------------------------------------------------------------------------
// Compute the Hull mass density summation.
//------------------------------------------------------------------------------

#include "computeHullSumMassDensity.hh"
#include "computeCSPHCorrections.hh"
#include "CSPHUtilities.hh"
#include "Field/FieldList.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Kernel/TableKernel.hh"
#include "NodeList/NodeList.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Utilities/comparisons.hh"

namespace Spheral {
namespace CSPHSpace {

using namespace std;
using std::min;
using std::max;
using std::abs;

using FieldSpace::FieldList;
using NeighborSpace::ConnectivityMap;
using KernelSpace::TableKernel;
using NodeSpace::NodeList;
using NodeSpace::FluidNodeList;
using BoundarySpace::Boundary;

namespace {
//------------------------------------------------------------------------------
// Worker methods specialized by dimension to compute the mass density based
// on convex hulls.
//------------------------------------------------------------------------------
// 1D
double hullMassDensity(const std::vector<Dim<1>::Vector>& pos,
                       const std::vector<Dim<1>::Scalar>& mass) {
  REQUIRE(pos.size() == mass.size());
  typedef Dim<1>::Scalar Scalar;
  typedef Dim<1>::Vector Vector;

  // Copy the two vectors to a single vector<pair> for sorting.
  vector<pair<Vector, Scalar> > stuff;
  const unsigned n = pos.size();
  for (unsigned i = 0; i != n; ++i) stuff.push_back(make_pair(pos[i], mass[i]));
  sort(stuff.begin(), stuff.end(), ComparePairByFirstElement<pair<Vector, Scalar> >());

  // Add up the masses.
  Scalar msum = 0.5*(stuff[0].second + stuff[n-1].second);
  for (unsigned i = 1; i < n - 1; ++i) msum += stuff[i].second;

  // Divide by volume and we're done.
  return msum*safeInv(stuff[n-1].first.x() - stuff[0].first.x());
}

// 2D
double hullMassDensity(const std::vector<Dim<2>::Vector>& pos,
                       const std::vector<Dim<2>::Scalar>& mass) {
  REQUIRE(pos.size() == mass.size());
  typedef Dim<2>::Scalar Scalar;
  typedef Dim<2>::Vector Vector;
  VERIFY(false);
}

// 3D
double hullMassDensity(const std::vector<Dim<3>::Vector>& pos,
                       const std::vector<Dim<3>::Scalar>& mass) {
  REQUIRE(pos.size() == mass.size());
  typedef Dim<3>::Scalar Scalar;
  typedef Dim<3>::Vector Vector;
  VERIFY(false);
}

}

//------------------------------------------------------------------------------
// The public method.
//------------------------------------------------------------------------------
template<typename Dimension>
void
computeHullSumMassDensity(const ConnectivityMap<Dimension>& connectivityMap,
                          const TableKernel<Dimension>& W,
                          const FieldList<Dimension, typename Dimension::Vector>& position,
                          const FieldList<Dimension, typename Dimension::Scalar>& mass,
                          const FieldList<Dimension, typename Dimension::SymTensor>& H,
                          FieldList<Dimension, typename Dimension::Scalar>& massDensity) {

  // Pre-conditions.
  const size_t numNodeLists = massDensity.size();
  REQUIRE(position.size() == numNodeLists);
  REQUIRE(mass.size() == numNodeLists);
  REQUIRE(H.size() == numNodeLists);

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename std::vector<BoundarySpace::Boundary<Dimension>*>::const_iterator ConstBoundaryIterator;

  const Scalar kernelExtent2 = FastMath::square(W.kernelExtent());

  // Walk the FluidNodeLists.
  for (size_t nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const FluidNodeList<Dimension>& nodeList = dynamic_cast<const FluidNodeList<Dimension>&>(massDensity[nodeListi]->nodeList());
    const int firstGhostNodei = nodeList.firstGhostNode();
    const Scalar rhoMin = nodeList.rhoMin();
    const Scalar rhoMax = nodeList.rhoMax();

    // Iterate over the nodes in this node list.
    for (typename ConnectivityMap<Dimension>::const_iterator iItr = connectivityMap.begin(nodeListi);
         iItr != connectivityMap.end(nodeListi);
         ++iItr) {
      const int i = *iItr;
      if (i < firstGhostNodei) {

        // Get the state for node i.
        const Vector& ri = position(nodeListi, i);
        const Scalar mi = mass(nodeListi, i);
        const SymTensor& Hi = H(nodeListi, i);

        // Get the neighbors for this node (in this NodeList).  We use the approximation here
        // that nodes from other NodeLists do not contribute to the density of this one.
        const vector<int>& connectivity = connectivityMap.connectivityForNode(nodeListi, i)[nodeListi];

        // Copy the neighbor positions & masses.
        vector<Vector> pos(1, ri);
        vector<Scalar> masses(1, mi);
        pos.reserve(connectivity.size() + 1);
        masses.reserve(connectivity.size() + 1);
        for (vector<int>::const_iterator jItr = connectivity.begin();
             jItr != connectivity.end();
             ++jItr) {
          const unsigned j = *jItr;
          const Vector& rj = position(nodeListi, j);
          const Scalar eta2i = (Hi*(ri - rj)).magnitude2();
          if (eta2i < kernelExtent2) {
            pos.push_back(position(nodeListi, *jItr));
            masses.push_back(mass(nodeListi, *jItr));
          }
        }
        CHECK(pos.size() >= 2);
        CHECK(masses.size() == pos.size());

        // Delegate to specialized methods.
        massDensity(nodeListi, i) = hullMassDensity(pos, masses);
        CHECK(massDensity(nodeListi, i) > 0.0);
      }
    }
  }
}

}
}

