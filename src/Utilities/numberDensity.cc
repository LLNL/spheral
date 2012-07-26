//---------------------------------Spheral++----------------------------------//
// numberDensity
//
// Compute the number density using the ASPH sum of the nodes.
//
// Created by JMO, Tue Feb  9 13:51:33 PST 2010
//----------------------------------------------------------------------------//
#include "numberDensity.hh"
#include "DataBase/DataBase.hh"
#include "Field/FieldList.hh"
#include "Kernel/TableKernel.hh"

namespace Spheral {

using namespace std;

using FieldSpace::FieldList;
using DataBaseSpace::DataBase;
using BoundarySpace::Boundary;
using KernelSpace::TableKernel;
using NeighborSpace::ConnectivityMap;
using NodeSpace::NodeList;

//------------------------------------------------------------------------------
// Return the number density for all nodes in the DataBase.
//------------------------------------------------------------------------------
template<typename Dimension>
FieldSpace::FieldList<Dimension, typename Dimension::Scalar>
numberDensity(const DataBaseSpace::DataBase<Dimension>& dataBase,
              const KernelSpace::TableKernel<Dimension>& W) {

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  // Get the state.
  const FieldList<Dimension, Vector> position = dataBase.globalPosition();
  const FieldList<Dimension, SymTensor> Hfield = dataBase.globalHfield();
  const ConnectivityMap<Dimension>& cm = dataBase.connectivityMap();
  const vector<const NodeList<Dimension>*>& nodeLists = cm.nodeLists();
  const size_t numNodeLists = dataBase.numNodeLists();

  // Prepare the result.
  FieldList<Dimension, Scalar> result = dataBase.newGlobalFieldList(0.0, "number density");

  // Iterate over the NodeLists.
  for (size_t nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {

    // Walk the nodes in this NodeList.
    for (typename ConnectivityMap<Dimension>::const_iterator iItr = cm.begin(nodeListi);
         iItr != cm.end(nodeListi); 
         ++iItr) {
      const int i = *iItr;
      const Vector& ri = position(nodeListi, i);
      const SymTensor& Hi = Hfield(nodeListi, i);
      const Scalar Hdeti = Hi.Determinant();
      const vector<vector<int> >& fullConnectivity = cm.connectivityForNode(nodeListi, i);
      CHECK(fullConnectivity.size() == numNodeLists);

      // Self-contribution.
      result(nodeListi, i) += W.kernelValue(0.0, Hdeti);

      // Walk the neighboring NodeLists.
      for (size_t nodeListj = 0; nodeListj != numNodeLists; ++nodeListj) {
        const int firstGhostNodej = nodeLists[nodeListj]->firstGhostNode();
        const vector<int>& connectivity = fullConnectivity[nodeListj];

        // Walk the neighbor in this NodeList.
        for (vector<int>::const_iterator jItr = connectivity.begin();
             jItr != connectivity.end();
             ++jItr) {
          const int j = *jItr;
          if (cm.calculatePairInteraction(nodeListi, i,
                                          nodeListj, j,
                                          firstGhostNodej)) {
            const Vector& rj = position(nodeListj, j);
            const SymTensor& Hj = Hfield(nodeListj, j);
            const Scalar Hdetj = Hj.Determinant();

            const Vector rij = ri - rj;
            const Scalar etai = (Hi*rij).magnitude();
            const Scalar etaj = (Hj*rij).magnitude();
            const Scalar Wi = W.kernelValue(etai, Hdeti);
            const Scalar Wj = W.kernelValue(etaj, Hdetj);

            result(nodeListi, i) += Wi;
            result(nodeListj, j) += Wj;
          }
        }
      }
    }
  }

  // That's it.
  return result;
}

}

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
namespace Spheral {

template FieldSpace::FieldList<Dim<1>, Dim<1>::Scalar> numberDensity<Dim<1> >(const DataBaseSpace::DataBase<Dim<1> >& dataBase, const KernelSpace::TableKernel<Dim<1> >& W);
template FieldSpace::FieldList<Dim<2>, Dim<2>::Scalar> numberDensity<Dim<2> >(const DataBaseSpace::DataBase<Dim<2> >& dataBase, const KernelSpace::TableKernel<Dim<2> >& W);
template FieldSpace::FieldList<Dim<3>, Dim<3>::Scalar> numberDensity<Dim<3> >(const DataBaseSpace::DataBase<Dim<3> >& dataBase, const KernelSpace::TableKernel<Dim<3> >& W);

}
