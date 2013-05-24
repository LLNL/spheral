//------------------------------------------------------------------------------
// Use SVPH to sample a FieldList.
//------------------------------------------------------------------------------
#include "sampleFieldListSVPH.hh"
#include "Field/FieldList.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Kernel/TableKernel.hh"
#include "NodeList/NodeList.hh"
#include "NodeList/FluidNodeList.hh"
#include "Utilities/safeInv.hh"

namespace Spheral {
namespace SVPHSpace {

using namespace std;
using std::min;
using std::max;
using std::abs;

using FieldSpace::FieldList;
using NeighborSpace::ConnectivityMap;
using KernelSpace::TableKernel;
using NodeSpace::NodeList;
using NodeSpace::FluidNodeList;

//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// The method itself.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
FieldList<Dimension, typename MathTraits<DataType>::GradientType>
sampleFieldListSVPH(const FieldList<Dimension, DataType>& field,
                    const FieldSpace::FieldList<Dimension, typename Dimension::Vector>& position,
                    const FieldSpace::FieldList<Dimension, typename Dimension::SymTensor>& H,
                    const NeighborSpace::ConnectivityMap<Dimension>& connectivityMap,
                    const KernelSpace::TableKernel<Dimension>& W,
                    const MeshSpace::Mesh<Dimension>& mesh,
                    const vector<Boundary<Dimension>*>& boundaries,
                    const bool firstOrderConsistent) {

  // Pre-conditions.
  const size_t numNodeLists = field.size();
  REQUIRE(position.size() == numNodeLists);
  REQUIRE(H.size() == numNodeLists);

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename MathTraits<DataType>::GradientType GradientType;

  // Prepare the result and some work fields.
  FieldList<Dimension, DataType> result(FieldList<Dimension, GradientType>::Copy);
  FieldList<Dimension, DataType> G(FieldList<Dimension, GradientType>::Copy);
  for (size_t nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const NodeList<Dimension>& nodeList = field[nodeListi]->nodeList();
    result.appendNewField("gradient of " + field[nodeListi]->name(), nodeList);
    G.appendNewField("internal gradient of " + field[nodeListi]->name(), nodeList);
  }

  // Make a first pass to evaluate the per cell gradients if we're doing first-order
  // consistency.
  if (firstOrderConsistent) {

    // Walk the FluidNodeLists.
    for (size_t nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
      const FluidNodeList<Dimension>& nodeList = dynamic_cast<const FluidNodeList<Dimension>&>(massDensity[nodeListi]->nodeList());
      const int firstGhostNodei = nodeList.firstGhostNode();

      

}

}
}
