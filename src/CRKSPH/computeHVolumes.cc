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
// Define a trait class to compute the H volume in each dimension.
//------------------------------------------------------------------------------
template<typename SymTensor> double Hvol(const SymTensor& H);

// 1D
template<> double Hvol<Dim<1>::SymTensor>(const Dim<1>::SymTensor& H) {
  return 2.0/H.xx();
}

// 2D
template<> double Hvol<Dim<2>::SymTensor>(const Dim<2>::SymTensor& H) {
  return M_PI/H.Determinant();
}

// 3D
template<> double Hvol<Dim<3>::SymTensor>(const Dim<3>::SymTensor& H) {
  return 4.0/3.0*M_PI/H.Determinant();
}

}

//------------------------------------------------------------------------------
// The volume method.
//------------------------------------------------------------------------------
template<typename Dimension>
void
computeHVolumes(const typename Dimension::Scalar kernelExtent,
                const FieldList<Dimension, typename Dimension::SymTensor>& H,
                FieldList<Dimension, typename Dimension::Scalar>& volume) {

  // Pre-conditions.
  const size_t numNodeLists = volume.size();
  REQUIRE(H.size() == numNodeLists);

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::SymTensor SymTensor;

  const Scalar kernelpow = Dimension::pownu(kernelExtent);

  // Walk the FluidNodeLists.
  for (size_t nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {

    // Iterate over the nodes in this node list.
    const unsigned n = volume[nodeListi]->numInternalElements();
    for (unsigned i = 0; i != n; ++i) {
      volume(nodeListi, i) = kernelpow * Hvol(H(nodeListi, i));
    }

  }
}

}

