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
// Return a MASH sampled version of the given FieldList at the new positions.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
FieldList<Dimension, DataType>
sampleFieldsMash(const FieldList<Dimension, DataType>& fieldList,
                 const FieldList<Dimension, typename Dimension::Vector>& position,
                 const FieldList<Dimension, typename Dimension::Scalar>& weight,
                 const FieldList<Dimension, typename Dimension::SymTensor>& Hfield,
                 const TableKernel<Dimension>& kernel,
                 const FieldList<Dimension, typename Dimension::Vector>& samplePositions,
                 const FieldList<Dimension, typename Dimension::Scalar>& sampleWeight,
                 const FieldList<Dimension, typename Dimension::SymTensor>& sampleHfield) {

  // Some convenient typedefs.
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  // Pre-conditions.
  CHECK(fieldList.numFields() == position.numFields());
  CHECK(fieldList.numFields() == weight.numFields());
  CHECK(fieldList.numFields() == Hfield.numFields());
  for (int i = 0; i < fieldList.numFields(); ++i) {
    CHECK(fieldList[i]->nodeListPtr() == position[i]->nodeListPtr());
    CHECK(fieldList[i]->nodeListPtr() == weight[i]->nodeListPtr());
    CHECK(fieldList[i]->nodeListPtr() == Hfield[i]->nodeListPtr());
  }
  CHECK(samplePositions.numFields() == sampleWeight.numFields());
  CHECK(samplePositions.numFields() == sampleHfield.numFields());
  for (int i = 0; i < samplePositions.numFields(); ++i) {
    CHECK(samplePositions[i]->nodeListPtr() == sampleWeight[i]->nodeListPtr());
    CHECK(samplePositions[i]->nodeListPtr() == sampleHfield[i]->nodeListPtr());
  }

  // Return FieldList.
  FieldList<Dimension, DataType> result;
  FieldList<Dimension, Scalar> normalization;
  result.copyFields();
  normalization.copyFields();
  for (typename FieldList<Dimension, Vector>::const_iterator
         fieldItr = samplePositions.begin();
       fieldItr < samplePositions.end(); 
       ++fieldItr) {
    result.appendField(Field<Dimension, DataType>("sample mash", (*fieldItr)->nodeList()));
    normalization.appendField(Field<Dimension, Scalar>("normalization", (*fieldItr)->nodeList(), 1.0e-30));
  }
  vector< vector<bool> > flagNodeDone(fieldList.numFields());
  for (typename FieldList<Dimension, DataType>::const_iterator
         fieldItr = fieldList.begin();
       fieldItr < fieldList.end(); 
       ++fieldItr) {
    flagNodeDone[fieldItr - fieldList.begin()].resize((*fieldItr)->nodeListPtr()->numInternalNodes(), false);
  }

  // Loop over all the positions in the fieldList we're sampling.
  for (InternalNodeIterator<Dimension> nodeItr = fieldList.internalNodeBegin();
       nodeItr < fieldList.internalNodeEnd();
       ++nodeItr) {

    // Check if this node has been done yet.
    if (!flagNodeDone[nodeItr.fieldID()][nodeItr.nodeID()]) {

      // Set the neighbor info over the positions we're sampling to.
      fieldList.setMasterNodeLists(position(nodeItr), Hfield(nodeItr));
      samplePositions.setMasterNodeLists(position(nodeItr), Hfield(nodeItr));

      // Loop over the set of master nodes in the FieldList we're sampling from.
      for (MasterNodeIterator<Dimension> masterItr = fieldList.masterNodeBegin();
           masterItr < fieldList.masterNodeEnd();
           ++masterItr) {
        CHECK(flagNodeDone[masterItr.fieldID()][masterItr.nodeID()] == false);
      
        // Sample node (i) state.
        const Vector& ri = position(masterItr);
        const SymTensor& Hi = Hfield(masterItr);
        const Scalar& weighti = weight(masterItr);
        const DataType& fieldi = fieldList(masterItr);

        // Refine the set of nodes we're sampling to for this position.
        samplePositions.setRefineNodeLists(ri, Hi);

        // Loop over the refined neighbors.
        for (RefineNodeIterator<Dimension> neighborItr = samplePositions.refineNodeBegin();
             neighborItr < samplePositions.refineNodeEnd();
             ++neighborItr) {

          // Node j's state.
          const Vector& rj = samplePositions(neighborItr);
          const SymTensor& Hj = sampleHfield(neighborItr);
          const Scalar& weightj = sampleWeight(neighborItr);

          const Vector rij = ri - rj;
          const Vector etai = Hi*rij;
          const Vector etaj = Hj*rij;
          CHECK(etai >= 0.0 && etaj >= 0.0);

          // Calculate the kernel estimates for each node.
          Scalar Wi = kernel(etai, 1.0);
          Scalar Wj = kernel(etaj, 1.0);

          // Get the symmetrized kernel weighting for this node pair.
          Scalar Wij, weightij;
          switch(neighborItr.nodeListPtr()->neighbor().neighborSearchType()) {
          case NeighborSpace::GatherScatter:
            Wij = 0.5*(Wi + Wj);
            weightij = 0.5*(weighti + weightj);
            break;

          case NeighborSpace::Gather:
            Wij = Wi;
            weightij = weighti;
            break;

          case NeighborSpace::Scatter:
            Wij = Wj;
            weightij = weightj;
            break;

          default:
            VERIFY2(false, "Unhandled neighbor search type.");
          }

          // Add this nodes contribution to the master value.
          normalization(neighborItr) += weightij*Wij;
          result(neighborItr) += fieldi*weightij*Wij;
        }
      
        // Flag this master node as done.
        flagNodeDone[masterItr.fieldID()][masterItr.nodeID()] = true;
      }
    }
  }

  // After we're done, all nodes in all NodeLists should be flagged as done.
  for (typename vector< vector<bool> >::const_iterator
         flagNodeItr = flagNodeDone.begin();
       flagNodeItr < flagNodeDone.end();
       ++flagNodeItr) {
    int checkcount = count(flagNodeItr->begin(), flagNodeItr->end(), false);
    if (checkcount > 0) {
      cerr << "Error in FieldList::sampleFieldsMash: Not all values determined on exit "
           << checkcount << endl;
    }
    CHECK(checkcount == 0);
  }

  // Normalize the estimate.
  result /= normalization;
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

//============================== sampleFieldsMash() ==============================
template 
FieldList<Dim<1>, Dim<1>::Scalar> 
sampleFieldsMash<Dim<1>, Dim<1>::Scalar>(const FieldList<Dim<1>, Dim<1>::Scalar>& fieldList,
                                         const FieldList<Dim<1>, Dim<1>::Vector>& position,
                                         const FieldList<Dim<1>, Dim<1>::Scalar>& weight,
                                         const FieldList<Dim<1>, Dim<1>::SymTensor>& Hfield,
                                         const TableKernel< Dim<1> >& kernel,
                                         const FieldList<Dim<1>, Dim<1>::Vector>& samplePositions,
                                         const FieldList<Dim<1>, Dim<1>::Scalar>& sampleWeight,
                                         const FieldList<Dim<1>, Dim<1>::SymTensor>& sampleHfield);

template 
FieldList<Dim<1>, Dim<1>::Vector> 
sampleFieldsMash<Dim<1>, Dim<1>::Vector>(const FieldList<Dim<1>, Dim<1>::Vector>& fieldList,
                                         const FieldList<Dim<1>, Dim<1>::Vector>& position,
                                         const FieldList<Dim<1>, Dim<1>::Scalar>& weight,
                                         const FieldList<Dim<1>, Dim<1>::SymTensor>& Hfield,
                                         const TableKernel< Dim<1> >& kernel,
                                         const FieldList<Dim<1>, Dim<1>::Vector>& samplePositions,
                                         const FieldList<Dim<1>, Dim<1>::Scalar>& sampleWeight,
                                         const FieldList<Dim<1>, Dim<1>::SymTensor>& sampleHfield);

template 
FieldList<Dim<1>, Dim<1>::Tensor> 
sampleFieldsMash<Dim<1>, Dim<1>::Tensor>(const FieldList<Dim<1>, Dim<1>::Tensor>& fieldList,
                                         const FieldList<Dim<1>, Dim<1>::Vector>& position,
                                         const FieldList<Dim<1>, Dim<1>::Scalar>& weight,
                                         const FieldList<Dim<1>, Dim<1>::SymTensor>& Hfield,
                                         const TableKernel< Dim<1> >& kernel,
                                         const FieldList<Dim<1>, Dim<1>::Vector>& samplePositions,
                                         const FieldList<Dim<1>, Dim<1>::Scalar>& sampleWeight,
                                         const FieldList<Dim<1>, Dim<1>::SymTensor>& sampleHfield);

template 
FieldList<Dim<1>, Dim<1>::SymTensor> 
sampleFieldsMash<Dim<1>, Dim<1>::SymTensor>(const FieldList<Dim<1>, Dim<1>::SymTensor>& fieldList,
                                            const FieldList<Dim<1>, Dim<1>::Vector>& position,
                                            const FieldList<Dim<1>, Dim<1>::Scalar>& weight,
                                            const FieldList<Dim<1>, Dim<1>::SymTensor>& Hfield,
                                            const TableKernel< Dim<1> >& kernel,
                                            const FieldList<Dim<1>, Dim<1>::Vector>& samplePositions,
                                            const FieldList<Dim<1>, Dim<1>::Scalar>& sampleWeight,
                                            const FieldList<Dim<1>, Dim<1>::SymTensor>& sampleHfield);

template 
FieldList<Dim<2>, Dim<2>::Scalar> 
sampleFieldsMash<Dim<2>, Dim<2>::Scalar>(const FieldList<Dim<2>, Dim<2>::Scalar>& fieldList,
                                         const FieldList<Dim<2>, Dim<2>::Vector>& position,
                                         const FieldList<Dim<2>, Dim<2>::Scalar>& weight,
                                         const FieldList<Dim<2>, Dim<2>::SymTensor>& Hfield,
                                         const TableKernel< Dim<2> >& kernel,
                                         const FieldList<Dim<2>, Dim<2>::Vector>& samplePositions,
                                         const FieldList<Dim<2>, Dim<2>::Scalar>& sampleWeight,
                                         const FieldList<Dim<2>, Dim<2>::SymTensor>& sampleHfield);

template 
FieldList<Dim<2>, Dim<2>::Vector> 
sampleFieldsMash<Dim<2>, Dim<2>::Vector>(const FieldList<Dim<2>, Dim<2>::Vector>& fieldList,
                                         const FieldList<Dim<2>, Dim<2>::Vector>& position,
                                         const FieldList<Dim<2>, Dim<2>::Scalar>& weight,
                                         const FieldList<Dim<2>, Dim<2>::SymTensor>& Hfield,
                                         const TableKernel< Dim<2> >& kernel,
                                         const FieldList<Dim<2>, Dim<2>::Vector>& samplePositions,
                                         const FieldList<Dim<2>, Dim<2>::Scalar>& sampleWeight,
                                         const FieldList<Dim<2>, Dim<2>::SymTensor>& sampleHfield);

template 
FieldList<Dim<2>, Dim<2>::Tensor> 
sampleFieldsMash<Dim<2>, Dim<2>::Tensor>(const FieldList<Dim<2>, Dim<2>::Tensor>& fieldList,
                                         const FieldList<Dim<2>, Dim<2>::Vector>& position,
                                         const FieldList<Dim<2>, Dim<2>::Scalar>& weight,
                                         const FieldList<Dim<2>, Dim<2>::SymTensor>& Hfield,
                                         const TableKernel< Dim<2> >& kernel,
                                         const FieldList<Dim<2>, Dim<2>::Vector>& samplePositions,
                                         const FieldList<Dim<2>, Dim<2>::Scalar>& sampleWeight,
                                         const FieldList<Dim<2>, Dim<2>::SymTensor>& sampleHfield);

template 
FieldList<Dim<2>, Dim<2>::SymTensor> 
sampleFieldsMash<Dim<2>, Dim<2>::SymTensor>(const FieldList<Dim<2>, Dim<2>::SymTensor>& fieldList,
                                            const FieldList<Dim<2>, Dim<2>::Vector>& position,
                                            const FieldList<Dim<2>, Dim<2>::Scalar>& weight,
                                            const FieldList<Dim<2>, Dim<2>::SymTensor>& Hfield,
                                            const TableKernel< Dim<2> >& kernel,
                                            const FieldList<Dim<2>, Dim<2>::Vector>& samplePositions,
                                            const FieldList<Dim<2>, Dim<2>::Scalar>& sampleWeight,
                                            const FieldList<Dim<2>, Dim<2>::SymTensor>& sampleHfield);

template 
FieldList<Dim<3>, Dim<3>::Scalar> 
sampleFieldsMash<Dim<3>, Dim<3>::Scalar>(const FieldList<Dim<3>, Dim<3>::Scalar>& fieldList,
                                         const FieldList<Dim<3>, Dim<3>::Vector>& position,
                                         const FieldList<Dim<3>, Dim<3>::Scalar>& weight,
                                         const FieldList<Dim<3>, Dim<3>::SymTensor>& Hfield,
                                         const TableKernel< Dim<3> >& kernel,
                                         const FieldList<Dim<3>, Dim<3>::Vector>& samplePositions,
                                         const FieldList<Dim<3>, Dim<3>::Scalar>& sampleWeight,
                                         const FieldList<Dim<3>, Dim<3>::SymTensor>& sampleHfield);

template 
FieldList<Dim<3>, Dim<3>::Vector> 
sampleFieldsMash<Dim<3>, Dim<3>::Vector>(const FieldList<Dim<3>, Dim<3>::Vector>& fieldList,
                                         const FieldList<Dim<3>, Dim<3>::Vector>& position,
                                         const FieldList<Dim<3>, Dim<3>::Scalar>& weight,
                                         const FieldList<Dim<3>, Dim<3>::SymTensor>& Hfield,
                                         const TableKernel< Dim<3> >& kernel,
                                         const FieldList<Dim<3>, Dim<3>::Vector>& samplePositions,
                                         const FieldList<Dim<3>, Dim<3>::Scalar>& sampleWeight,
                                         const FieldList<Dim<3>, Dim<3>::SymTensor>& sampleHfield);

template 
FieldList<Dim<3>, Dim<3>::Tensor> 
sampleFieldsMash<Dim<3>, Dim<3>::Tensor>(const FieldList<Dim<3>, Dim<3>::Tensor>& fieldList,
                                         const FieldList<Dim<3>, Dim<3>::Vector>& position,
                                         const FieldList<Dim<3>, Dim<3>::Scalar>& weight,
                                         const FieldList<Dim<3>, Dim<3>::SymTensor>& Hfield,
                                         const TableKernel< Dim<3> >& kernel,
                                         const FieldList<Dim<3>, Dim<3>::Vector>& samplePositions,
                                         const FieldList<Dim<3>, Dim<3>::Scalar>& sampleWeight,
                                         const FieldList<Dim<3>, Dim<3>::SymTensor>& sampleHfield);

template 
FieldList<Dim<3>, Dim<3>::SymTensor> 
sampleFieldsMash<Dim<3>, Dim<3>::SymTensor>(const FieldList<Dim<3>, Dim<3>::SymTensor>& fieldList,
                                            const FieldList<Dim<3>, Dim<3>::Vector>& position,
                                            const FieldList<Dim<3>, Dim<3>::Scalar>& weight,
                                            const FieldList<Dim<3>, Dim<3>::SymTensor>& Hfield,
                                            const TableKernel< Dim<3> >& kernel,
                                            const FieldList<Dim<3>, Dim<3>::Vector>& samplePositions,
                                            const FieldList<Dim<3>, Dim<3>::Scalar>& sampleWeight,
                                            const FieldList<Dim<3>, Dim<3>::SymTensor>& sampleHfield);

}
}
