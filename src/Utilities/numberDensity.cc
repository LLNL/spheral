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

using std::vector;
using std::string;
using std::pair;
using std::make_pair;
using std::cout;
using std::cerr;
using std::endl;
using std::min;
using std::max;
using std::abs;

namespace Spheral {

//------------------------------------------------------------------------------
// Return the number density for all nodes in the DataBase.
//------------------------------------------------------------------------------
template<typename Dimension>
FieldList<Dimension, typename Dimension::Scalar>
numberDensity(const DataBase<Dimension>& dataBase,
              const TableKernel<Dimension>& W) {

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

