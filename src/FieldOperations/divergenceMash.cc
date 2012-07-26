//---------------------------------Spheral++----------------------------------//
// FieldListFunctionsMash -- A set of global functions which can be applied to
// FieldLists using MASH prescriptions.
//
// Created by JMO, Wed Dec  6 21:09:29 PST 2000
//----------------------------------------------------------------------------//

#include "FieldListFunctionsMash.hh"
#include "Field/FieldList.hh"
#include "Field/Field.hh"
#include "Field/NodeIterators.hh"
#include "NodeList/NodeList.hh"
#include "Neighbor/Neighbor.hh"
#include "Kernel/TableKernel.hh"
#include "Geometry/MathTraits.hh"

namespace Spheral {
namespace FieldSpace {

using namespace std;
using NodeSpace::NodeList;
using NeighborSpace::Neighbor;
using KernelSpace::TableKernel;

//------------------------------------------------------------------------------
// Calculate the divergence of a FieldList.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
FieldList<Dimension, typename MathTraits<Dimension, DataType>::DivergenceType>
divergenceMash(const FieldList<Dimension, DataType>& fieldList,
               const FieldList<Dimension, typename Dimension::Vector>& position,
               const FieldList<Dimension, typename Dimension::Scalar>& weight,
               const FieldList<Dimension, typename Dimension::SymTensor>& Hfield,
               const TableKernel<Dimension>& kernel) {

  // Typedef's to ease typing/understandability.
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename MathTraits<Dimension, DataType>::DivergenceType DivergenceType;

  // Return FieldList.
  FieldList<Dimension, DivergenceType> result;
  vector< vector<bool> > flagNodeDone(fieldList.numFields());
  result.copyFields();
  for (typename FieldList<Dimension, DataType>::const_iterator fieldItr = fieldList.begin();
       fieldItr < fieldList.end(); 
       ++fieldItr) {
    result.appendField(Field<Dimension, DivergenceType>("div mash", (*fieldItr)->nodeList()));
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

        // Loop over the refined neighbors.
        const Vector& ri = position(masterItr);
        const SymTensor& Hi = Hfield(masterItr);
        const Scalar& weighti = weight(masterItr);
        const DataType& fieldi = fieldList(masterItr);

        Scalar normalization = 0.0;
        for (RefineNodeIterator<Dimension> neighborItr = fieldList.refineNodeBegin();
             neighborItr < fieldList.refineNodeEnd();
             ++neighborItr) {

          const Vector& rj = position(neighborItr);
          const SymTensor& Hj = Hfield(neighborItr);
          const Scalar& weightj = weight(neighborItr);
          const DataType& fieldj = fieldList(neighborItr);

          const Vector rij = ri - rj;
          const Vector etai = Hi*rij;
          const Vector etaj = Hj*rij;
          const Vector etaiNorm = etai.unitVector();
          const Vector etajNorm = etaj.unitVector();
          CHECK(etai >= 0.0 && etaj >= 0.0);

          // Get the symmetrized kernel gradient for this node pair.
          Scalar Wij;
          Vector gradWij;
          switch((*fieldList.begin())->nodeListPtr()->neighbor().neighborSearchType()) {
          case NeighborSpace::GatherScatter:
            Wij = 0.5*(kernel(etai, 1.0) + kernel(etaj, 1.0));
            gradWij = 0.5*(Hi*etaiNorm*kernel.grad(etai, 1.0) + 
                           Hj*etajNorm*kernel.grad(etaj, 1.0));
            break;

          case NeighborSpace::Gather:
            Wij = kernel(etai, 1.0);
            gradWij = Hi*etaiNorm*kernel.grad(etai, 1.0);
            break;

          case NeighborSpace::Scatter:
            Wij = kernel(etaj, 1.0);
            gradWij = Hj*etajNorm*kernel.grad(etaj, 1.0);
            break;

          default:
            VERIFY2(false, "Unhandled neighbor search type.");
          }

          // Add this nodes contribution to the master value.
          normalization += weightj*Wij;
          result(masterItr) += weightj*(fieldj - fieldi).dot(gradWij);

          //           Scalar Wi = kernel(etai, Hi);
          //           Scalar Wj = kernel(etaj, Hj);
          //           Scalar Wij = 0.5*(Wi + Wj);
          //           result(masterItr) += 2.0*
          //             Wj*((weighti*fieldi - weightj*fieldj)/(rij.magnitude2() + FLT_MIN)).dot(rij);
        }

        // Normalize the divergence value.
        CHECK(normalization > 0.0);
        result(masterItr) /= normalization;

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
      cerr << "Error in FieldList::smoothFieldsMash: Not all values determined on exit "
           << checkcount << endl;
    }
    CHECK(checkcount == 0);
  }

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

//============================== divergenceMash() ==============================
template 
FieldList<Dim<1>, MathTraits<Dim<1>, Dim<1>::Vector>::DivergenceType> 
divergenceMash<Dim<1>, Dim<1>::Vector>(const FieldList<Dim<1>, Dim<1>::Vector>& fieldList,
                                       const FieldList<Dim<1>, Dim<1>::Vector>& position,
                                       const FieldList<Dim<1>, Dim<1>::Scalar>& weight,
                                       const FieldList<Dim<1>, Dim<1>::SymTensor>& Hfield,
                                       const TableKernel< Dim<1> >& kernel);

template 
FieldList<Dim<1>, MathTraits<Dim<1>, Dim<1>::Tensor>::DivergenceType> 
divergenceMash<Dim<1>, Dim<1>::Tensor>(const FieldList<Dim<1>, Dim<1>::Tensor>& fieldList,
                                       const FieldList<Dim<1>, Dim<1>::Vector>& position,
                                       const FieldList<Dim<1>, Dim<1>::Scalar>& weight,
                                       const FieldList<Dim<1>, Dim<1>::SymTensor>& Hfield,
                                       const TableKernel< Dim<1> >& kernel);
template 
FieldList<Dim<2>, MathTraits<Dim<2>, Dim<2>::Vector>::DivergenceType> 
divergenceMash<Dim<2>, Dim<2>::Vector>(const FieldList<Dim<2>, Dim<2>::Vector>& fieldList,
                                       const FieldList<Dim<2>, Dim<2>::Vector>& position,
                                       const FieldList<Dim<2>, Dim<2>::Scalar>& weight,
                                       const FieldList<Dim<2>, Dim<2>::SymTensor>& Hfield,
                                       const TableKernel< Dim<2> >& kernel);

template 
FieldList<Dim<2>, MathTraits<Dim<2>, Dim<2>::Tensor>::DivergenceType> 
divergenceMash<Dim<2>, Dim<2>::Tensor>(const FieldList<Dim<2>, Dim<2>::Tensor>& fieldList,
                                       const FieldList<Dim<2>, Dim<2>::Vector>& position,
                                       const FieldList<Dim<2>, Dim<2>::Scalar>& weight,
                                       const FieldList<Dim<2>, Dim<2>::SymTensor>& Hfield,
                                       const TableKernel< Dim<2> >& kernel);
template 
FieldList<Dim<3>, MathTraits<Dim<3>, Dim<3>::Vector>::DivergenceType> 
divergenceMash<Dim<3>, Dim<3>::Vector>(const FieldList<Dim<3>, Dim<3>::Vector>& fieldList,
                                       const FieldList<Dim<3>, Dim<3>::Vector>& position,
                                       const FieldList<Dim<3>, Dim<3>::Scalar>& weight,
                                       const FieldList<Dim<3>, Dim<3>::SymTensor>& Hfield,
                                       const TableKernel< Dim<3> >& kernel);

template 
FieldList<Dim<3>, MathTraits<Dim<3>, Dim<3>::Tensor>::DivergenceType> 
divergenceMash<Dim<3>, Dim<3>::Tensor>(const FieldList<Dim<3>, Dim<3>::Tensor>& fieldList,
                                       const FieldList<Dim<3>, Dim<3>::Vector>& position,
                                       const FieldList<Dim<3>, Dim<3>::Scalar>& weight,
                                       const FieldList<Dim<3>, Dim<3>::SymTensor>& Hfield,
                                       const TableKernel< Dim<3> >& kernel);
}
}
