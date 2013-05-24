//------------------------------------------------------------------------------
// Use SVPH to sample a FieldList.
//------------------------------------------------------------------------------
#include "sampleFieldListSVPH.hh"
#include "NodeList/NodeList.hh"
#include "NodeList/FluidNodeList.hh"
#include "Utilities/safeInv.hh"
#include "Geometry/Dimension.hh"
#include "Geometry/MathTraits.hh"

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
using BoundarySpace::Boundary;

namespace {
//------------------------------------------------------------------------------
// Convert 1D types to scalars.
//------------------------------------------------------------------------------
inline Dim<1>::Scalar scalarAdapter(const Dim<1>::Scalar& val) { return val; }
inline Dim<1>::Scalar scalarAdapter(const Dim<1>::Vector& val) { return val.x(); }
inline Dim<1>::Scalar scalarAdapter(const Dim<1>::Tensor& val) { return val.xx(); }
inline Dim<1>::Scalar scalarAdapter(const Dim<1>::SymTensor& val) { return val.xx(); }
  
}

//------------------------------------------------------------------------------
// 1-D gradient in cell.
//------------------------------------------------------------------------------
template<typename DataType>
typename MathTraits<Dim<1>, DataType>::GradientType
computeCellGradient(const FieldSpace::FieldList<Dim<1>, DataType>& fieldList,
                    const FieldSpace::FieldList<Dim<1>, Dim<1>::Vector>& position,
                    const int nodeListi,
                    const int i,
                    const MeshSpace::Mesh<Dim<1> >& mesh) {
  using namespace FieldSpace;
  using namespace MeshSpace;
  typedef Dim<1> Dimension;
  typedef Dimension::Scalar Scalar;
  typedef Dimension::Vector Vector;
  typedef Mesh<Dimension>::Zone Zone;
  typedef Mesh<Dimension>::Face Face;
  typedef typename MathTraits<Dim<1>, DataType>::GradientType GradientType;

  // Find the neighboring mesh cells.
  const Zone& zonei = mesh.zone(nodeListi, i);
  const std::vector<int>& faces = zonei.faceIDs();
  CHECK(faces.size() == 2);
  const int face1ID = Mesh<Dimension>::positiveID(faces[0]),
            face2ID = Mesh<Dimension>::positiveID(faces[1]);
  const Face face1 = mesh.face(face1ID),
             face2 = mesh.face(face2ID);
  const int zone1ID = Mesh<Dimension>::positiveID(face1.oppositeZoneID(zonei.ID())),
            zone2ID = Mesh<Dimension>::positiveID(face2.oppositeZoneID(zonei.ID()));

  // Look up the indices into the FieldList for the two neighboring cells.
  unsigned nodeList1ID = nodeListi, node1ID = i, 
           nodeList2ID = nodeListi, node2ID = i;
  if (zone1ID != Mesh<Dimension>::UNSETID) mesh.lookupNodeListID(zone1ID, nodeList1ID, node1ID);
  if (zone2ID != Mesh<Dimension>::UNSETID) mesh.lookupNodeListID(zone2ID, nodeList2ID, node2ID);

  // Now compute the weighted intrazonal gradient.
  const Scalar fi = scalarAdapter(fieldList(nodeListi, i)),
               f1 = scalarAdapter(fieldList(nodeList1ID, node1ID)),
               f2 = scalarAdapter(fieldList(nodeList2ID, node2ID));
  const Scalar dx01 = (position(nodeList1ID, node1ID) - position(nodeListi, i)).x(),
               dx02 = (position(nodeList2ID, node2ID) - position(nodeListi, i)).x();
  CHECK(std::abs(dx01) + std::abs(dx02) > 0.0);
  const GradientType result((sgn(dx01)*(f1 - fi) + sgn(dx02)*(f2 - fi))*
                            safeInv(std::abs(dx01) + std::abs(dx02)));
  return result;
}

//------------------------------------------------------------------------------
// The method itself.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
FieldList<Dimension, DataType>
sampleFieldListSVPH(const FieldList<Dimension, DataType>& fieldList,
                    const FieldSpace::FieldList<Dimension, typename Dimension::Vector>& position,
                    const FieldSpace::FieldList<Dimension, typename Dimension::SymTensor>& Hfield,
                    const NeighborSpace::ConnectivityMap<Dimension>& connectivityMap,
                    const KernelSpace::TableKernel<Dimension>& W,
                    const MeshSpace::Mesh<Dimension>& mesh,
                    const vector<Boundary<Dimension>*>& boundaries,
                    const bool firstOrderConsistent) {

  // Pre-conditions.
  const size_t numNodeLists = fieldList.size();
  REQUIRE(position.size() == numNodeLists);
  REQUIRE(Hfield.size() == numNodeLists);

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename MathTraits<Dimension, DataType>::GradientType GradientType;

  // Prepare the result and some work fields.
  FieldList<Dimension, DataType> result(FieldList<Dimension, DataType>::Copy);
  FieldList<Dimension, GradientType> G(FieldList<Dimension, GradientType>::Copy);
  for (size_t nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const NodeList<Dimension>& nodeList = fieldList[nodeListi]->nodeList();
    result.appendNewField("SVPH sample of " + fieldList[nodeListi]->name(), nodeList, DataType(0));
    G.appendNewField("internal gradient of " + fieldList[nodeListi]->name(), nodeList, GradientType(0));
  }

  // Make a first pass to evaluate the per cell gradients if we're doing first-order
  // consistency.
  if (firstOrderConsistent) {

    // Walk the NodeLists.
    for (int nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
      const NodeList<Dimension>& nodeList = fieldList[nodeListi]->nodeList();

      // Iterate over the nodes in this node list.
      for (int i = 0; i != nodeList.numInternalNodes(); ++i) {

        // Do the dimension specific lookup of the internal cell gradient, based on the cell geometry.
        G(nodeListi, i) = computeCellGradient(fieldList, position, nodeListi, i, mesh);
      }
    }
  }

  // Walk the NodeLists to build our answer.
  const Scalar W0 = W.kernelValue(0.0, 1.0);
  for (int nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const NodeList<Dimension>& nodeList = fieldList[nodeListi]->nodeList();

    // Iterate over the nodes in this node list.
    for (typename ConnectivityMap<Dimension>::const_iterator iItr = connectivityMap.begin(nodeListi);
         iItr != connectivityMap.end(nodeListi);
         ++iItr) {
      const int i = *iItr;

      // Get the state for node i.
      const Vector& ri = position(nodeListi, i);
      const Scalar Vi = mesh.zone(nodeListi, i).volume();
      Scalar norm = Vi*W0;
      result(nodeListi, i) = Vi*W0 * fieldList(nodeListi, i);

      // Walk the neighbors for this node.
      const vector<vector<int> >& fullConnectivity = connectivityMap.connectivityForNode(nodeListi, i);
      for (int nodeListj = 0; nodeListj != numNodeLists; ++nodeListj) {
        const vector<int>& connectivity = fullConnectivity[nodeListj];
        for (vector<int>::const_iterator jItr = connectivity.begin();
             jItr != connectivity.end();
             ++jItr) {
          const int j = *jItr;

          // Get the state for node j.
          const DataType& Fj = fieldList(nodeListj, j);
          const Vector& rj = position(nodeListj, j);
          const SymTensor& Hj = Hfield(nodeListj, j);
          const Scalar Vj = mesh.zone(nodeListj, j).volume();
          const GradientType& Gj = G(nodeListj, j);

          // Pair-wise kernel type stuff.
          const Vector rij = ri - rj;
          const Vector etaj = Hj*rij;
          const Scalar Wj = W.kernelValue(etaj.magnitude(), 1.0);

          // Increment the result.
          norm += Vj*Wj;
          result(nodeListi, i) += Vj*Wj * (Fj + Gj.dot(rij));
        }
      }

      // Normalize our ith value.
      CHECK(norm > 0.0);
      result(nodeListi, i) /= norm;
    }
  }

  // That's it.
  return result;
}

}
}
