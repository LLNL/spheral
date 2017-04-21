//---------------------------------Spheral++------------------------------------
// Compute a volume per point based on the local H tensor.
//------------------------------------------------------------------------------
#include "computeHVolumes.hh"
#include "Field/Field.hh"
#include "Field/FieldList.hh"
#include "NodeList/NodeList.hh"

namespace Spheral {

using namespace std;

using FieldSpace::Field;
using FieldSpace::FieldList;
using NodeSpace::NodeList;

namespace {
//------------------------------------------------------------------------------
// Define a dimension specialized method to compute the H volume.
//------------------------------------------------------------------------------
// 1D
inline
double Hvol(const Dim<1>::SymTensor& H,
            const Dim<1>::Scalar nPerh) {
  return 1.0/(nPerh*H.xx());
}

// 2D
inline
double Hvol(const Dim<2>::SymTensor& H,
            const Dim<2>::Scalar nPerh) {
  return 0.25*M_PI/(nPerh*H).Determinant();
}

// 3D
inline
double Hvol(const Dim<3>::SymTensor& H,
            const Dim<3>::Scalar nPerh) {
  return 0.125*4.0/3.0*M_PI/(nPerh*H).Determinant();
}

}

//------------------------------------------------------------------------------
// The volume method.
//------------------------------------------------------------------------------
template<typename Dimension>
void
computeHVolumes(const typename Dimension::Scalar nPerh,
                const FieldList<Dimension, typename Dimension::SymTensor>& H,
                FieldList<Dimension, typename Dimension::Scalar>& volume) {

  // Pre-conditions.
  const size_t numNodeLists = volume.size();
  REQUIRE(nPerh > 0.0);
  REQUIRE(H.size() == numNodeLists);

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::SymTensor SymTensor;

  // Walk the FluidNodeLists.
  for (size_t nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {

    // Iterate over the nodes in this node list.
    const unsigned n = volume[nodeListi]->numInternalElements();
    for (unsigned i = 0; i != n; ++i) {
      volume(nodeListi, i) = Hvol(H(nodeListi, i), nPerh);
    }
  }
}

}

