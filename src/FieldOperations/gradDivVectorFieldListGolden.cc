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
#include "Field/FieldList.hh"
#include "Field/Field.hh"
#include "Field/NodeIterators.hh"
#include "NodeList/NodeList.hh"
#include "Neighbor/Neighbor.hh"
#include "Kernel/TableKernel.hh"
#include "Boundary/Boundary.hh"

#include "TAU.h"

namespace Spheral {
namespace FieldSpace {

using namespace std;
using NodeSpace::NodeList;
using NeighborSpace::Neighbor;
using KernelSpace::TableKernel;
using BoundarySpace::Boundary;

//------------------------------------------------------------------------------
// Calculate the gradient of the divergence of a Vector FieldList.
// Golden rule method.
//------------------------------------------------------------------------------
template<typename Dimension>
FieldList<Dimension, typename Dimension::Vector>
gradDivVectorFieldListGolden
(const FieldList<Dimension, typename Dimension::Vector>& fieldList,
 const FieldList<Dimension, typename Dimension::Vector>& position,
 const FieldList<Dimension, typename Dimension::Scalar>& weight,
 const FieldList<Dimension, typename Dimension::Scalar>& mass,
 const FieldList<Dimension, typename Dimension::Scalar>& rho,
 const FieldList<Dimension, typename Dimension::SymTensor>& Hfield,
 const TableKernel<Dimension>& kernel) {

  // TAU timers.
  TAU_PROFILE("gradDivVectorFieldListGolden", "", TAU_USER);
  TAU_PROFILE_TIMER(TimeGoldenSecondReserve, "gradDivVectorFieldListGolden", "Reserve return Fields", TAU_USER);
  TAU_PROFILE_TIMER(TimeGoldenSecondMaster, "gradDivVectorFieldListGolden", "Select master nodes", TAU_USER);
  TAU_PROFILE_TIMER(TimeGoldenSecondRefine, "gradDivVectorFieldListGolden", "Select refine nodes", TAU_USER);
  TAU_PROFILE_TIMER(TimeGoldenSecondNodeI, "gradDivVectorFieldListGolden", "State for node I", TAU_USER);
  TAU_PROFILE_TIMER(TimeGoldenSecondNodeJ, "gradDivVectorFieldListGolden", "State for node J", TAU_USER);
  TAU_PROFILE_TIMER(TimeGoldenSecondRefineLoop, "gradDivVectorFieldListGolden", "Main loop over refine neighobors", TAU_USER);
  TAU_PROFILE_TIMER(TimeGoldenSecondW, "gradDivVectorFieldListGolden", "Calculate W, gradW, grad2W", TAU_USER);
  TAU_PROFILE_TIMER(TimeGoldenSecondAccumulate, "gradDivVectorFieldListGolden", "Accumulate components", TAU_USER);
  TAU_PROFILE_TIMER(TimeGoldenSecondFinalize, "gradDivVectorFieldListGolden", "Determine final answer", TAU_USER);
  TAU_PROFILE_TIMER(TimeGoldenSecondCheck, "gradDivVectorFieldListGolden", "Completeness check", TAU_USER);

  // Some convenient typedefs.
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  // Return FieldList.
  TAU_PROFILE_START(TimeGoldenSecondReserve);
  FieldList<Dimension, Vector> result;
  vector< vector<bool> > flagNodeDone(fieldList.numFields());
  result.copyFields();
  for (typename FieldList<Dimension, Vector>::const_iterator fieldItr = fieldList.begin();
       fieldItr < fieldList.end(); 
       ++fieldItr) {
    result.appendField(Field<Dimension, Vector>("grad div " + (*fieldItr)->name(), (*fieldItr)->nodeList()));
    flagNodeDone[fieldItr - fieldList.begin()].resize((*fieldItr)->nodeListPtr()->numInternalNodes(), false);
  }
  TAU_PROFILE_STOP(TimeGoldenSecondReserve);

  // Remember the square of the kernel extent.
  const Scalar cutoff2 = kernel.kernelExtent()*kernel.kernelExtent();

  // Loop over all the elements in the input FieldList.
  for (InternalNodeIterator<Dimension> nodeItr = fieldList.internalNodeBegin();
       nodeItr < fieldList.internalNodeEnd();
       ++nodeItr) {

    // Check if this node has been done yet.
    if (!flagNodeDone[nodeItr.fieldID()][nodeItr.nodeID()]) {

      // We will do the batch of master nodes associated with this node together.
      // Set the neighbor information.
      TAU_PROFILE_START(TimeGoldenSecondMaster);
      fieldList.setMasterNodeLists(position(nodeItr), Hfield(nodeItr));
      TAU_PROFILE_STOP(TimeGoldenSecondMaster);

      // Now loop over all the master nodes.
      for (MasterNodeIterator<Dimension> masterItr = fieldList.masterNodeBegin();
           masterItr < fieldList.masterNodeEnd();
           ++masterItr) {
        CHECK(flagNodeDone[masterItr.fieldID()][masterItr.nodeID()] == false);

        // Set the refined neighbor information for this master node.
        TAU_PROFILE_START(TimeGoldenSecondRefine);
        fieldList.setRefineNodeLists(position(masterItr), Hfield(masterItr));
        TAU_PROFILE_STOP(TimeGoldenSecondRefine);

        // State for this node.
        TAU_PROFILE_START(TimeGoldenSecondNodeI);
        const Vector& ri = position(masterItr);
        const SymTensor& Hi = Hfield(masterItr);
        const Scalar& mi = mass(masterItr);
        const Scalar& rhoi = rho(masterItr);
        const Scalar& weighti = weight(masterItr);
        const Vector& fieldi = fieldList(masterItr);
        CHECK(Hi.Determinant() > 0.0);
        CHECK(mi > 0.0);
        CHECK(rhoi > 0.0);
        TAU_PROFILE_STOP(TimeGoldenSecondNodeI);

        // Temp variables to accumulate the grad field and grad rho for this
        // node.
        vector<Tensor> grad2elements(Dimension::nDim);
        Vector gradRho;
        Tensor gradF;

        // Loop over the refined neighbors, and calculate the various
        // sums that contributed to grad div F.
        for (RefineNodeIterator<Dimension> neighborItr = fieldList.refineNodeBegin();
             neighborItr < fieldList.refineNodeEnd();
             ++neighborItr) {
          TAU_PROFILE_START(TimeGoldenSecondRefineLoop);
          if (neighborItr != masterItr) {

//             TAU_PROFILE_START((TimeGoldenSecondNodeJ);
            const Vector& rj = position(neighborItr);
            const SymTensor& Hj = Hfield(neighborItr);
            const Scalar& mj = mass(neighborItr);
            const Scalar& rhoj = rho(neighborItr);
            const Scalar& weightj = weight(neighborItr);
            const Vector& fieldj = fieldList(neighborItr);
            CHECK(Hj.Determinant() > 0.0);
            CHECK(mj > 0.0);
            CHECK(rhoj > 0.0);

            const Vector rij = ri - rj;
            const Vector etai = Hi*rij;
            const Vector etaj = Hj*rij;
//             TAU_PROFILE_STOP(TimeGoldenSecondNodeJ);

            if (etai.magnitude2() < cutoff2 ||
                etaj.magnitude2() < cutoff2) {
//               TAU_PROFILE_START(TimeGoldenSecondW);
              const Vector etaiNorm = etai.unitVector();
              const Vector etajNorm = etaj.unitVector();
              const Vector Hetai = Hi*etaiNorm;
              const Vector Hetaj = Hj*etajNorm;

              const Scalar getai = kernel.grad(etai, Hi);
              const Scalar getaj = kernel.grad(etaj, Hj);
              const Scalar g2etai = kernel.grad2(etai, Hi);
              const Scalar g2etaj = kernel.grad2(etaj, Hj);

              const Vector gWi = Hetai*getai;
              const Vector gWj = Hetaj*getaj;

              const Tensor acki = Hetai.dyad(Hetai);
              const Tensor ackj = Hetaj.dyad(Hetaj);

              const Tensor H2i = Hi*Hi;
              const Tensor H2j = Hj*Hj;

              const Tensor g2Wi = acki*g2etai + 
                (H2i.Transpose() - acki)/(etai.magnitude() + 1.0e-30)*getai;
              const Tensor g2Wj = ackj*g2etaj +
                (H2j.Transpose() - ackj)/(etaj.magnitude() + 1.0e-30)*getaj;

              // Get the symmetrized kernel weighting for this node pair.
              Vector gWij;
              Tensor g2Wij;
              switch((*fieldList.begin())->nodeListPtr()->neighbor().neighborSearchType()) {
              case NeighborSpace::GatherScatter:
                gWij = 0.5*(gWi + gWj);
                g2Wij = 0.5*(g2Wi + g2Wj);
                break;

              case NeighborSpace::Gather:
                gWij = gWi;
                g2Wij = g2Wi;
                break;

              case NeighborSpace::Scatter:
                gWij = gWj;
                g2Wij = g2Wj;
                break;
              }
//               TAU_PROFILE_STOP(TimeGoldenSecondW);

              // Sum the contributions.
//               TAU_PROFILE_START(TimeGoldenSecondAccumulate);

//               const Scalar weightij = 0.5*(weighti + weightj);
//               const Vector fji = weightij*(fieldj - fieldi);
//               for (int gamma = 0; gamma < Dimension::nDim; ++gamma) {
//                 for (int beta = 0; beta < Dimension::nDim; ++beta) {
//                   for (int alpha = 0; alpha < Dimension::nDim; ++alpha) {
//                     grad2elements[gamma](beta, alpha) += fji(gamma)*g2Wij(beta, alpha);
//                   }
//                   gradF(gamma, beta) += fji(gamma)*gWij(beta);
//                 }
//                 gradRho(gamma) += weightij*(rhoj - rhoi)*gWij(gamma);
//               }

              const Vector fij = fieldi - fieldj;
              for (int gamma = 0; gamma < Dimension::nDim; ++gamma) {
                for (int beta = 0; beta < Dimension::nDim; ++beta) {
                  for (int alpha = 0; alpha < Dimension::nDim; ++alpha) {
                    grad2elements[gamma](beta, alpha) -= mj*fij(gamma)*g2Wij(beta, alpha);
                  }
                  gradF(gamma, beta) -= mj*fij(gamma)*gWij(beta);
                }
                gradRho(gamma) -= mj*(rhoi - rhoj)*gWij(gamma);
              }
//               TAU_PROFILE_STOP(TimeGoldenSecondAccumulate);
            }
          }
          TAU_PROFILE_STOP(TimeGoldenSecondRefineLoop);
        }

        // Put together the complete derivative.
        TAU_PROFILE_START(TimeGoldenSecondFinalize);
        CHECK(rhoi > 0.0);
        gradRho /= rhoi;
        gradF /= rhoi;
        for (int beta = 0; beta < Dimension::nDim; ++beta) {
          for (int alpha = 0; alpha < Dimension::nDim; ++alpha) {
            result(masterItr)(beta) += grad2elements[alpha](beta, alpha) -
              gradRho(alpha)*gradF(alpha, beta) -
              gradRho(beta)*gradF(alpha, alpha);
          }
        }
        result(masterItr) /= rhoi;
        TAU_PROFILE_STOP(TimeGoldenSecondFinalize);

        // This master node is finished.
        flagNodeDone[masterItr.fieldID()][masterItr.nodeID()] = true;
      }
    }
  }

  // After we're done, all nodes in all NodeLists should be flagged as done.
  TAU_PROFILE_START(TimeGoldenSecondCheck);
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
  TAU_PROFILE_STOP(TimeGoldenSecondCheck);

  return result;
}

}
}

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
namespace Spheral {
namespace FieldSpace {

using KernelSpace::TableKernel;

//========================== gradDivVectorFieldListGolden() ==========================
template 
FieldList<Dim<1>, Dim<1>::Vector> 
gradDivVectorFieldListGolden< Dim<1> >
(const FieldList<Dim<1>, Dim<1>::Vector>& fieldList,
 const FieldList<Dim<1>, Dim<1>::Vector>& position,
 const FieldList<Dim<1>, Dim<1>::Scalar>& weight,
 const FieldList<Dim<1>, Dim<1>::Scalar>& mass,
 const FieldList<Dim<1>, Dim<1>::Scalar>& rho,
 const FieldList<Dim<1>, Dim<1>::SymTensor>& Hfield,
 const TableKernel< Dim<1> >& kernel);

template 
FieldList<Dim<2>, Dim<2>::Vector> 
gradDivVectorFieldListGolden< Dim<2> >
(const FieldList<Dim<2>, Dim<2>::Vector>& fieldList,
 const FieldList<Dim<2>, Dim<2>::Vector>& position,
 const FieldList<Dim<2>, Dim<2>::Scalar>& weight,
 const FieldList<Dim<2>, Dim<2>::Scalar>& mass,
 const FieldList<Dim<2>, Dim<2>::Scalar>& rho,
 const FieldList<Dim<2>, Dim<2>::SymTensor>& Hfield,
 const TableKernel< Dim<2> >& kernel);

template 
FieldList<Dim<3>, Dim<3>::Vector> 
gradDivVectorFieldListGolden< Dim<3> >
(const FieldList<Dim<3>, Dim<3>::Vector>& fieldList,
 const FieldList<Dim<3>, Dim<3>::Vector>& position,
 const FieldList<Dim<3>, Dim<3>::Scalar>& weight,
 const FieldList<Dim<3>, Dim<3>::Scalar>& mass,
 const FieldList<Dim<3>, Dim<3>::Scalar>& rho,
 const FieldList<Dim<3>, Dim<3>::SymTensor>& Hfield,
 const TableKernel< Dim<3> >& kernel);

}
}
