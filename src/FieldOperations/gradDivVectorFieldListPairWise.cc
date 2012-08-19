//---------------------------------Spheral++----------------------------------//
// FieldListSecondDerivatives
// A set of experimental methods for evaluating the second derivative of a 
// FieldList.
//
// Created by JMO, Wed Dec 18 22:46:54 PST 2002
//----------------------------------------------------------------------------//

#include <vector>
using std::vector;

#include "FieldListSecondDerivatives.hh"
#include "PairWiseFieldListFunctions.hh"
#include "Field/FieldList.hh"
#include "Field/Field.hh"
#include "Field/NodeIterators.hh"
#include "NodeList/NodeList.hh"
#include "Neighbor/Neighbor.hh"
#include "Kernel/TableKernel.hh"
#include "Boundary/Boundary.hh"
#include "Utilities/rotationMatrix.hh"

namespace Spheral {
namespace FieldSpace {

using namespace std;
using NodeSpace::NodeList;
using NeighborSpace::Neighbor;
using KernelSpace::TableKernel;
using BoundarySpace::Boundary;

//------------------------------------------------------------------------------
// Calculate the gradient of the divergence of a Vector FieldList.
// This version uses a direct pairwise method to evaluate the first derivative,
// and applies a kernel gradient in the same sum to evaluate the second 
// derivative.
//------------------------------------------------------------------------------
template<typename Dimension>
FieldList<Dimension, typename Dimension::Vector>
gradDivVectorFieldListPairWise
(const FieldList<Dimension, typename Dimension::Vector>& fieldList,
 const FieldList<Dimension, typename Dimension::Vector>& position,
 const FieldList<Dimension, typename Dimension::Scalar>& weight,
 const FieldList<Dimension, typename Dimension::Scalar>& mass,
 const FieldList<Dimension, typename Dimension::Scalar>& rho,
 const FieldList<Dimension, typename Dimension::SymTensor>& Hfield,
 const TableKernel<Dimension>& kernel) {

  // Some convenient typedefs.
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  // Return FieldList.
  FieldList<Dimension, Vector> result;
  vector< vector<bool> > flagNodeDone(fieldList.numFields());
  result.copyFields();
  for (typename FieldList<Dimension, Vector>::const_iterator fieldItr = fieldList.begin();
       fieldItr < fieldList.end(); 
       ++fieldItr) {
    result.appendField(Field<Dimension, Vector>("grad div " + (*fieldItr)->name(), (*fieldItr)->nodeList()));
    flagNodeDone[fieldItr - fieldList.begin()].resize((*fieldItr)->nodeListPtr()->numInternalNodes(), false);
  }

  // Loop over all the elements in the input FieldList.
  for (InternalNodeIterator<Dimension> nodeItr = fieldList.internalNodeBegin();
       nodeItr < fieldList.internalNodeEnd();
       ++nodeItr) {

    // Check if this node has been done yet.
    if (!flagNodeDone[nodeItr.fieldID()][nodeItr.nodeID()]) {

      // We will do the batch of master nodes associated with this node together.
      // Set the neighbor information.
      fieldList.setMasterNodeLists(position(nodeItr), Hfield(nodeItr));

      // Now loop over all the master nodes.
      for (MasterNodeIterator<Dimension> masterItr = fieldList.masterNodeBegin();
           masterItr < fieldList.masterNodeEnd();
           ++masterItr) {
        CHECK(flagNodeDone[masterItr.fieldID()][masterItr.nodeID()] == false);

        // Set the refined neighbor information for this master node.
        fieldList.setRefineNodeLists(position(masterItr), Hfield(masterItr));

        // State for node i.
        const Vector& ri = position(masterItr);
        const SymTensor& Hi = Hfield(masterItr);
        const Scalar& weighti = weight(masterItr);
        const Scalar& mi = mass(masterItr);
        const Scalar& rhoi = rho(masterItr);
        const Vector& fieldi = fieldList(masterItr);

        // Loop over the refined neighbors.
        vector<Tensor> grad2i(Dimension::nDim);
        Tensor normalization(0.0);
        for (RefineNodeIterator<Dimension> neighborItr = fieldList.refineNodeBegin();
             neighborItr < fieldList.refineNodeEnd();
             ++neighborItr) {
          if (neighborItr != masterItr) {

            const Vector& rj = position(neighborItr);
            const SymTensor& Hj = Hfield(neighborItr);
            const Scalar& weightj = weight(neighborItr);
            const Scalar& mj = mass(neighborItr);
            const Scalar& rhoj = rho(neighborItr);
            const Vector& fieldj = fieldList(neighborItr);

            const Vector rij = ri - rj;
            const Vector etai = Hi*rij;
            const Vector etaj = Hj*rij;

            const Vector etaiNorm = etai.unitVector();
            const Vector etajNorm = etaj.unitVector();
            const Vector Hetai = Hi*etaiNorm;
            const Vector Hetaj = Hj*etajNorm;

            const Scalar Wi = kernel(etai, Hi);
            const Scalar Wj = kernel(etaj, Hj);

            const Scalar getai = kernel.grad(etai, Hi);
            const Scalar getaj = kernel.grad(etaj, Hj);

            const Vector gWi = Hetai*getai;
            const Vector gWj = Hetaj*getaj;

            // Get the symmetrized kernel weighting for this node pair.
            Scalar Wij;
            Vector gWij;
            switch((*fieldList.begin())->nodeListPtr()->neighbor().neighborSearchType()) {
            case NeighborSpace::GatherScatter:
              Wij = 0.5*(Wi + Wj);
              gWij = 0.5*(gWi + gWj);
              break;

            case NeighborSpace::Gather:
              Wij = Wi;
              gWij = gWi;
              break;

            case NeighborSpace::Scatter:
              Wij = Wj;
              gWij = gWj;
              break;

            default:
              VERIFY2(false, "Unhandled neighbor search type.");
            }

            // Sum this pairs contribution to the elements.
            const Vector fij = fieldi - fieldj;
            const Scalar wij = weightj*Wij;
            const Vector gwij = weightj*gWij;

            const Scalar frach = (0.01/Hi.Trace())/Dimension::nDim;
            const Vector rjiUnit = -rij.unitVector();
            const Tensor R = rotationMatrix(rjiUnit);
            const Tensor Rinverse = R.Transpose();
            const Scalar dxp = rij.magnitude()/(rij.magnitude2() + frach*frach);

            // Represent this pairs delta as the first column of a matrix
            // in the rotated frame.
            const Vector dfdxpcol = -dxp*R*fij;
            Tensor dfdxp;
            dfdxp.setColumn(0, dfdxpcol);

            // Rotate the delta to the lab frame.
            dfdxp.rotationalTransform(Rinverse);

            // Increment the result with the delta.
            for (int i = 0; i < Dimension::nDim; ++i) {
              grad2i[i] += dfdxp*gwij(i);
            }

            // Update the normalization.
            Tensor tweight(0.0);
            for (int i = 0; i < Dimension::nDim; ++i) tweight(i,0) = wij;
            tweight.rotationalTransform(Rinverse);
            for (int i = 0; i < Dimension::nDim; ++i) {
              for (int j = 0; j < Dimension::nDim; ++j) {
                normalization(i,j) += std::abs(tweight(i,j));
              }
            }

//             // Sum this pairs contribution to the elements.
//             const Vector fij = fieldi - fieldj;
//             const Scalar divfij = fij.dot(rij)/(rij.magnitude2() + 1e-4);
//             normalization += weightj*Wij;
//             result(masterItr) += weightj*divfij*gWij;

          }
        }

        // Apply the normalization factor.
//         CHECK(normalization > 0.0);
//         result(masterItr) *= Dimension::nDim/normalization;
        for (int i = 0; i < Dimension::nDim; ++i) {
          for (int j = 0; j < Dimension::nDim; ++j) {
            for (int k = 0; k < Dimension::nDim; ++k) {
              grad2i[i](j,k) /= normalization(j,k) + 1.0e-10;
            }
          }
        }

        // Now contract of the rank 3 tensor to get the divergence.
        for (int i = 0; i < Dimension::nDim; ++i) {
          result(masterItr)(i) = grad2i[i].Trace();
        }
        
        // This master node is finished.
        flagNodeDone[masterItr.fieldID()][masterItr.nodeID()] = true;

      }
    }
  }

  // After we're done, all nodes in all NodeLists should be flagged as done.
  for (typename vector< vector<bool> >::const_iterator flagNodeItr = flagNodeDone.begin();
       flagNodeItr < flagNodeDone.end();
       ++flagNodeItr) {
    int checkcount = count(flagNodeItr->begin(), flagNodeItr->end(), false);
    if (checkcount > 0) {
      cerr << "Error in FieldList::smoothFields: Not all values determined on exit "
           << checkcount << endl;
    }
    CHECK(checkcount == 0);
  }

  return result;
}



// //------------------------------------------------------------------------------
// // Calculate the gradient of the divergence of a Vector FieldList.
// // Explicit method that performs the div and grad operations as sequential first
// // derivatives (in this case using the pair wise operators).
// //------------------------------------------------------------------------------
// template<typename Dimension>
// FieldList<Dimension, typename Dimension::Vector>
// gradDivVectorFieldListPairWise
// (const FieldList<Dimension, typename Dimension::Vector>& fieldList,
//  const FieldList<Dimension, typename Dimension::Vector>& position,
//  const FieldList<Dimension, typename Dimension::Scalar>& weight,
//  const FieldList<Dimension, typename Dimension::Scalar>& mass,
//  const FieldList<Dimension, typename Dimension::Scalar>& density,
//  const FieldList<Dimension, typename Dimension::SymTensor>& Hfield,
//  const TableKernel<Dimension>& kernel,
//  const vector<Boundary<Dimension>*>&  boundaries) {

//   // Some convenient typedefs.
//   typedef typename Dimension::Scalar Scalar;
//   typedef typename Dimension::Vector Vector;
//   typedef typename Dimension::Tensor Tensor;
//   typedef typename Dimension::SymTensor SymTensor;

//   // First evaluate the divergence of the input field list.
//   FieldList<Dimension, Scalar> divField = divergencePairWise(fieldList,
//                                                              position,
//                                                              weight,
//                                                              mass,
//                                                              density,
//                                                              Hfield,
//                                                              kernel);

//   // Apply boundary conditions to the divergence.
//   for (typename vector<Boundary<Dimension>*>::const_iterator bcItr = boundaries.begin();
//        bcItr < boundaries.end();
//        ++bcItr) {
//     (*bcItr)->applyFieldListGhostBoundary(divField);
//   }

//   // Now take the gradient of this.
//   FieldList<Dimension, Vector> result = gradientPairWise(divField,
//                                                          position,
//                                                          weight,
//                                                          mass,
//                                                          density,
//                                                          Hfield,
//                                                          kernel);

//   // That's it.
//   return result;
// }

}
}

