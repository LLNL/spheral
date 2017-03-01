//---------------------------------Spheral++----------------------------------//
// FieldListSecondDerivatives
// A set of experimental methods for evaluating the second derivative of a 
// FieldList.
//
// Created by JMO, Wed Dec 18 22:46:54 PST 2002
//----------------------------------------------------------------------------//

#include <vector>

#include "FieldListSecondDerivatives.hh"
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
// Simplest method, just use the second derivative of the kernel.
//------------------------------------------------------------------------------
template<typename Dimension>
FieldList<Dimension, typename Dimension::Vector>
gradDivVectorFieldListSimple
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
        Scalar Wnormalization = weighti*kernel(0.0, Hi);
        vector<Tensor> D2fDx2i(Dimension::nDim);
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

            const Vector rji = rj - ri;
            const Vector etai = Hi*rji;
            const Vector etaj = Hj*rji;

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
            case NeighborSpace::NeighborSearchType::GatherScatter:
              Wij = 0.5*(Wi + Wj);
              gWij = 0.5*(gWi + gWj);
              break;

            case NeighborSpace::NeighborSearchType::Gather:
              Wij = Wi;
              gWij = gWi;
              break;

            case NeighborSpace::NeighborSearchType::Scatter:
              Wij = Wj;
              gWij = gWj;
              break;

            default:
              VERIFY2(false, "Unhandled neighbor search type.");
            }

            // Sum this pairs contribution to the elements.
            Wnormalization += weightj*Wij;
            const Vector fji = fieldj - fieldi;
            const Scalar frach = 0.01/Hi.Trace()/Dimension::nDim;
            const Vector rjiUnit = rji.unitVector();
            const Tensor R = rotationMatrix(rjiUnit);
            const Tensor Rinverse = R.Transpose();
            const Scalar dxp = rji.magnitude()/(rji.magnitude2() + frach*frach);

            // Find the local contribution to dfdx.
            const Vector dfdxcol = (R*fji)*dxp;

            // Take the second derivative in the prime frame.
            const Tensor d2fdx2 = dfdxcol*(R*gWij);

            // Now increment the total result in the prime frame, and then rotate
            // it back.
            result(masterItr) = R*result(masterItr);
            result(masterItr) -= weightj*(d2fdx2.diagonalElements());
            result(masterItr) = Rinverse*result(masterItr);
          }
        }

        // Put together the elements into the final grad div answer.
        CHECK(Wnormalization > 0.0);
        result(masterItr) *= 2.0/Wnormalization;

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

// template<typename Dimension>
// FieldList<Dimension, typename Dimension::Vector>
// gradDivVectorFieldListSimple
// (const FieldList<Dimension, typename Dimension::Vector>& fieldList,
//  const FieldList<Dimension, typename Dimension::Vector>& position,
//  const FieldList<Dimension, typename Dimension::Scalar>& weight,
//  const FieldList<Dimension, typename Dimension::Scalar>& mass,
//  const FieldList<Dimension, typename Dimension::Scalar>& rho,
//  const FieldList<Dimension, typename Dimension::SymTensor>& Hfield,
//  const TableKernel<Dimension>& kernel) {

//   // Some convenient typedefs.
//   typedef typename Dimension::Scalar Scalar;
//   typedef typename Dimension::Vector Vector;
//   typedef typename Dimension::Tensor Tensor;
//   typedef typename Dimension::SymTensor SymTensor;

//   // Return FieldList.
//   FieldList<Dimension, Vector> result;
//   vector< vector<bool> > flagNodeDone(fieldList.numFields());
//   result.copyFields();
//   for (typename FieldList<Dimension, Vector>::const_iterator fieldItr = fieldList.begin();
//        fieldItr < fieldList.end(); 
//        ++fieldItr) {
//     result.appendField(Field<Dimension, Vector>((*fieldItr)->nodeList()));
//     flagNodeDone[fieldItr - fieldList.begin()].resize((*fieldItr)->nodeListPtr()->numInternalNodes(), false);
//   }

//   // Loop over all the elements in the input FieldList.
//   for (InternalNodeIterator<Dimension> nodeItr = fieldList.internalNodeBegin();
//        nodeItr < fieldList.internalNodeEnd();
//        ++nodeItr) {

//     // Check if this node has been done yet.
//     if (!flagNodeDone[nodeItr.fieldID()][nodeItr.nodeID()]) {

//       // We will do the batch of master nodes associated with this node together.
//       // Set the neighbor information.
//       fieldList.setMasterNodeLists(position(nodeItr), Hfield(nodeItr));

//       // Now loop over all the master nodes.
//       for (MasterNodeIterator<Dimension> masterItr = fieldList.masterNodeBegin();
//            masterItr < fieldList.masterNodeEnd();
//            ++masterItr) {
//         CHECK(flagNodeDone[masterItr.fieldID()][masterItr.nodeID()] == false);

//         // Set the refined neighbor information for this master node.
//         fieldList.setRefineNodeLists(position(masterItr), Hfield(masterItr));

//         // State for node i.
//         const Vector& ri = position(masterItr);
//         const SymTensor& Hi = Hfield(masterItr);
//         const Scalar& mi = mass(masterItr);
//         const Scalar& rhoi = rho(masterItr);
//         const Vector& fieldi = fieldList(masterItr);

//         // Loop over the refined neighbors again, and calculate grad div f.
//         vector<Tensor> grad2elements(Dimension::nDim);
//         for (RefineNodeIterator<Dimension> neighborItr = fieldList.refineNodeBegin();
//              neighborItr < fieldList.refineNodeEnd();
//              ++neighborItr) {
//           if (neighborItr != masterItr) {

//             const Vector& rj = position(neighborItr);
//             const SymTensor& Hj = Hfield(neighborItr);
//             const Scalar& weightj = weight(neighborItr);
//             const Scalar& mj = mass(neighborItr);
//             const Scalar& rhoj = rho(neighborItr);
//             const Vector& fieldj = fieldList(neighborItr);

//             const Vector rji = rj - ri;
//             const Vector etai = Hi*rji;
//             const Vector etaj = Hj*rji;

//             const Vector etaiNorm = etai.unitVector();
//             const Vector etajNorm = etaj.unitVector();
//             const Vector Hetai = Hi*etaiNorm;
//             const Vector Hetaj = Hj*etajNorm;

//             const Scalar getai = kernel.grad(etai, Hi);
//             const Scalar getaj = kernel.grad(etaj, Hj);
//             const Scalar g2etai = kernel.grad2(etai, Hi);
//             const Scalar g2etaj = kernel.grad2(etaj, Hj);

//             const Vector gWi = Hetai*getai;
//             const Vector gWj = Hetaj*getaj;

//             const Tensor acki = Hetai.dyad(Hetai);
//             const Tensor ackj = Hetaj.dyad(Hetaj);

//             const Tensor H2i = Hi*Hi;
//             const Tensor H2j = Hj*Hj;

//             const Tensor g2Wi = acki*g2etai +
//               (H2i.Transpose() - acki)/(etai.magnitude() + 1.0e-10)*getai;
//             const Tensor g2Wj = ackj*g2etaj +
//               (H2j.Transpose() - ackj)/(etaj.magnitude() + 1.0e-10)*getaj;

//             // Get the symmetrized kernel weighting for this node pair.
//             Vector gWij;
//             Tensor g2Wij;
//             switch((*fieldList.begin())->nodeListPtr()->neighbor().neighborSearchType()) {
//             case Neighbor<Dimension>::GatherScatter:
//               gWij = 0.5*(gWi + gWj);
//               g2Wij = 0.5*(g2Wi + g2Wj);
//               break;

//             case Neighbor<Dimension>::Gather:
//               gWij = gWi;
//               g2Wij = g2Wi;
//               break;

//             case Neighbor<Dimension>::Scatter:
//               gWij = gWj;
//               g2Wij = g2Wj;
//               break;
//             }

//             // Sum this pairs contribution to the elements.
//             const Vector fji = fieldj - fieldi;
//             for (int gamma = 0; gamma < Dimension::nDim; ++gamma) {
//               for (int beta = 0; beta < Dimension::nDim; ++beta) {
//                 for (int alpha = 0; alpha < Dimension::nDim; ++alpha) {
//                   grad2elements[gamma](beta, alpha) +=
//                     weightj*weightj*fieldj(gamma)*g2Wij(beta, alpha);
//                 }
//               }
//             }
//           }
//         }

//         // Put together the elements into the final grad div answer.
//         for (int beta = 0; beta < Dimension::nDim; ++beta) {
//           for (int alpha = 0; alpha < Dimension::nDim; ++alpha) {
//             result(masterItr)(beta) += grad2elements[alpha](beta, alpha);
//           }
//         }

//         // This master node is finished.
//         flagNodeDone[masterItr.fieldID()][masterItr.nodeID()] = true;

//       }
//     }
//   }

//   // After we're done, all nodes in all NodeLists should be flagged as done.
//   for (typename vector< vector<bool> >::const_iterator flagNodeItr = flagNodeDone.begin();
//        flagNodeItr < flagNodeDone.end();
//        ++flagNodeItr) {
//     int checkcount = count(flagNodeItr->begin(), flagNodeItr->end(), false);
//     if (checkcount > 0) {
//       cerr << "Error in FieldList::smoothFields: Not all values determined on exit "
//            << checkcount << endl;
//     }
//     CHECK(checkcount == 0);
//   }

//   return result;
// }

}
}

