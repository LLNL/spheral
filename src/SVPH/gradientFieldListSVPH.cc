//------------------------------------------------------------------------------
// Use SVPH to take the gradient of a FieldList.
//------------------------------------------------------------------------------
#include "sampleFieldListSVPH.hh"
#include "computeSVPHCorrections.hh"
#include "NodeList/NodeList.hh"
#include "NodeList/FluidNodeList.hh"
#include "Utilities/safeInv.hh"
#include "Geometry/Dimension.hh"
#include "Geometry/innerProduct.hh"

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
using Geometry::innerProduct;

//------------------------------------------------------------------------------
// Local utility methods.
//------------------------------------------------------------------------------
namespace {

// 1D
Dim<1>::Vector
graddot(const Dim<1>::Vector& xij, const Dim<1>::Vector& Bi, const Dim<1>::Tensor& gradBi) {
  return Dim<1>::Vector(Bi(0) + xij(0)*gradBi(0,0));
}

// 2D
Dim<2>::Vector
graddot(const Dim<2>::Vector& xij, const Dim<2>::Vector& Bi, const Dim<2>::Tensor& gradBi) {
  return Dim<2>::Vector(Bi(0) + xij(0)*gradBi(0,0) + xij(1)*gradBi(0,1) + xij(1)*(gradBi(1,0) - gradBi(0,1)),
                        Bi(1) + xij(0)*gradBi(1,0) + xij(1)*gradBi(1,1) - xij(0)*(gradBi(1,0) - gradBi(0,1)));
}

// 3D
Dim<3>::Vector
graddot(const Dim<3>::Vector& xij, const Dim<3>::Vector& Bi, const Dim<3>::Tensor& gradBi) {
  return Dim<3>::Vector();
}

}

//------------------------------------------------------------------------------
// The method itself.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
FieldList<Dimension, typename MathTraits<Dimension, DataType>::GradientType>
gradientFieldListSVPH(const FieldList<Dimension, DataType>& fieldList,
                      const FieldSpace::FieldList<Dimension, typename Dimension::Vector>& position,
                      const FieldSpace::FieldList<Dimension, typename Dimension::SymTensor>& Hfield,
                      const NeighborSpace::ConnectivityMap<Dimension>& connectivityMap,
                      const KernelSpace::TableKernel<Dimension>& W,
                      const MeshSpace::Mesh<Dimension>& mesh,
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
  FieldList<Dimension, GradientType> result(FieldList<Dimension, GradientType>::Copy);
  FieldList<Dimension, Scalar> volume(FieldList<Dimension, Scalar>::Copy);
  FieldList<Dimension, Vector> B(FieldList<Dimension, Vector>::Copy);
  FieldList<Dimension, Tensor> gradB(FieldList<Dimension, Tensor>::Copy);
  for (size_t nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const NodeList<Dimension>& nodeList = fieldList[nodeListi]->nodeList();
    result.appendNewField("SVPH gradient of " + fieldList[nodeListi]->name(), nodeList, GradientType());
    volume.appendNewField("volume", nodeList, 0.0);
    B.appendNewField("SVPH linear correction for " + fieldList[nodeListi]->name(), nodeList, Vector::zero);
    gradB.appendNewField("SVPH linear correction gradient for " + fieldList[nodeListi]->name(), nodeList, Tensor::zero);
  }

  // If we're enforcing first-order consistency, compute the correction fields.
  if (firstOrderConsistent) {

    // Copy the cell volumes to the FieldList.
    for (int nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
      const NodeList<Dimension>& nodeList = fieldList[nodeListi]->nodeList();
      for (int i = 0; i != nodeList.numNodes(); ++i) {
        volume(nodeListi, i) = mesh.zone(nodeListi, i).volume();
      }
    }

    // We have a handy utility method to compute the corrections.
    computeSVPHCorrections(connectivityMap,
                           W,
                           volume,
                           position,
                           Hfield,
                           B,
                           gradB);
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
      const Vector& Bi = B(nodeListi, i);
      const Tensor& gradBi = gradB(nodeListi, i);
      const DataType& Fi = fieldList(nodeListi, i);
      Scalar norm = Vi*W0;

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

          // Pair-wise kernel type stuff.
          const Vector rij = ri - rj;
          const Vector etaj = Hj*rij;
          const Vector Hetaj = Hj*etaj.unitVector();
          const pair<double, double> WWj = W.kernelAndGradValue(etaj.magnitude(), 1.0);
          const Scalar Wj = WWj.first;
          const Scalar gWj = WWj.second;
          const Vector gradWj = gWj*Hetaj;

          // Increment the result.
          norm += Vj*(1.0 + Bi.dot(rij))*Wj;
          result(nodeListi, i) += Vj*(Fj - Fi)*((1.0 + Bi.dot(rij))*gradWj + graddot(rij, Bi, gradBi)*Wj);
          // result(nodeListi, i) += Vj*(Fj - Fi)*((1.0 + Bi.dot(rij))*gradWj + (Bi + innerProduct<Dimension>(rij, gradBi))*Wj);
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
