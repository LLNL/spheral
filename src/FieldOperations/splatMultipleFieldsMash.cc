//---------------------------------Spheral++----------------------------------//
// FieldListFunctionsMash -- A set of global functions which can be applied to
// FieldLists using MASH prescriptions.
//
// Created by JMO, Tue Jun  4 17:07:31 PDT 2002
//----------------------------------------------------------------------------//
#include "FieldListFunctionsMash.hh"
#include "Field/FieldList.hh"
#include "Field/FieldListSet.hh"
#include "Field/Field.hh"
#include "Field/NodeIterators.hh"
#include "NodeList/NodeList.hh"
#include "Neighbor/Neighbor.hh"
#include "Kernel/TableKernel.hh"
#include "Geometry/MathTraits.hh"

#include "DBC.hh"

namespace Spheral {
namespace FieldSpace {

using namespace std;
using NodeSpace::NodeList;
using NeighborSpace::Neighbor;
using KernelSpace::TableKernel;

//------------------------------------------------------------------------------
// Return a MASH donated version of the given FieldList at the new positions.
// This version uses the "splat" operation, guaranteeing that this operation is
// conservative in the sense that summing over the original Field values and
// the new Field values yeilds the same total.
// This version operates on a set of FieldLists at the same time, purely as an 
// optimzation over the single FieldList version, splatFieldsMash.
//------------------------------------------------------------------------------
template<typename Dimension>
FieldListSet<Dimension>
splatMultipleFieldsMash(const FieldListSet<Dimension>& fieldListSet,
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
  for (typename vector< FieldList<Dimension, Scalar> >::const_iterator fieldListItr = 
         fieldListSet.ScalarFieldLists.begin();
       fieldListItr < fieldListSet.ScalarFieldLists.end();
       ++fieldListItr) {
    const FieldList<Dimension, Scalar>& fieldList = *fieldListItr;
    REQUIRE(fieldList.numFields() == position.numFields());
    REQUIRE(fieldList.numFields() == weight.numFields());
    REQUIRE(fieldList.numFields() == Hfield.numFields());
    for (int i = 0; i < fieldList.numFields(); ++i) {
      REQUIRE(fieldList[i]->nodeListPtr() == position[i]->nodeListPtr());
      REQUIRE(fieldList[i]->nodeListPtr() == weight[i]->nodeListPtr());
      REQUIRE(fieldList[i]->nodeListPtr() == Hfield[i]->nodeListPtr());
    }
  }
  for (typename vector< FieldList<Dimension, Vector> >::const_iterator fieldListItr = 
         fieldListSet.VectorFieldLists.begin();
       fieldListItr < fieldListSet.VectorFieldLists.end();
       ++fieldListItr) {
    const FieldList<Dimension, Vector>& fieldList = *fieldListItr;
    REQUIRE(fieldList.numFields() == position.numFields());
    REQUIRE(fieldList.numFields() == weight.numFields());
    REQUIRE(fieldList.numFields() == Hfield.numFields());
    for (int i = 0; i < fieldList.numFields(); ++i) {
      REQUIRE(fieldList[i]->nodeListPtr() == position[i]->nodeListPtr());
      REQUIRE(fieldList[i]->nodeListPtr() == weight[i]->nodeListPtr());
      REQUIRE(fieldList[i]->nodeListPtr() == Hfield[i]->nodeListPtr());
    }
  }
  for (typename vector< FieldList<Dimension, Tensor> >::const_iterator fieldListItr = 
         fieldListSet.TensorFieldLists.begin();
       fieldListItr < fieldListSet.TensorFieldLists.end();
       ++fieldListItr) {
    const FieldList<Dimension, Tensor>& fieldList = *fieldListItr;
    REQUIRE(fieldList.numFields() == position.numFields());
    REQUIRE(fieldList.numFields() == weight.numFields());
    REQUIRE(fieldList.numFields() == Hfield.numFields());
    for (int i = 0; i < fieldList.numFields(); ++i) {
      REQUIRE(fieldList[i]->nodeListPtr() == position[i]->nodeListPtr());
      REQUIRE(fieldList[i]->nodeListPtr() == weight[i]->nodeListPtr());
      REQUIRE(fieldList[i]->nodeListPtr() == Hfield[i]->nodeListPtr());
    }
  }
  for (typename vector< FieldList<Dimension, SymTensor> >::const_iterator fieldListItr = 
         fieldListSet.SymTensorFieldLists.begin();
       fieldListItr < fieldListSet.SymTensorFieldLists.end();
       ++fieldListItr) {
    const FieldList<Dimension, SymTensor>& fieldList = *fieldListItr;
    REQUIRE(fieldList.numFields() == position.numFields());
    REQUIRE(fieldList.numFields() == weight.numFields());
    REQUIRE(fieldList.numFields() == Hfield.numFields());
    for (int i = 0; i < fieldList.numFields(); ++i) {
      REQUIRE(fieldList[i]->nodeListPtr() == position[i]->nodeListPtr());
      REQUIRE(fieldList[i]->nodeListPtr() == weight[i]->nodeListPtr());
      REQUIRE(fieldList[i]->nodeListPtr() == Hfield[i]->nodeListPtr());
    }
  }
  REQUIRE(samplePositions.numFields() == sampleWeight.numFields());
  REQUIRE(samplePositions.numFields() == sampleHfield.numFields());
  for (int i = 0; i < samplePositions.numFields(); ++i) {
    REQUIRE(samplePositions[i]->nodeListPtr() == sampleWeight[i]->nodeListPtr());
    REQUIRE(samplePositions[i]->nodeListPtr() == sampleHfield[i]->nodeListPtr());
  }

  // Return FieldList.
  FieldListSet<Dimension> resultSet;
  resultSet.ScalarFieldLists.resize(fieldListSet.ScalarFieldLists.size());
  resultSet.VectorFieldLists.resize(fieldListSet.VectorFieldLists.size());
  resultSet.TensorFieldLists.resize(fieldListSet.TensorFieldLists.size());
  resultSet.SymTensorFieldLists.resize(fieldListSet.SymTensorFieldLists.size());
  for (typename vector<FieldList<Dimension, Scalar> >::iterator resultItr = 
         resultSet.ScalarFieldLists.begin();
       resultItr < resultSet.ScalarFieldLists.end();
       ++resultItr) resultItr->copyFields();
  for (typename vector<FieldList<Dimension, Vector> >::iterator resultItr = 
         resultSet.VectorFieldLists.begin();
       resultItr < resultSet.VectorFieldLists.end();
       ++resultItr) resultItr->copyFields();
  for (typename vector<FieldList<Dimension, Tensor> >::iterator resultItr = 
         resultSet.TensorFieldLists.begin();
       resultItr < resultSet.TensorFieldLists.end();
       ++resultItr) resultItr->copyFields();
  for (typename vector<FieldList<Dimension, SymTensor> >::iterator resultItr = 
         resultSet.SymTensorFieldLists.begin();
       resultItr < resultSet.SymTensorFieldLists.end();
       ++resultItr) resultItr->copyFields();
  for (typename FieldList<Dimension, Vector>::const_iterator fieldItr = samplePositions.begin();
       fieldItr < samplePositions.end(); 
       ++fieldItr) {
    for (typename vector<FieldList<Dimension, Scalar> >::iterator resultItr = 
           resultSet.ScalarFieldLists.begin();
         resultItr < resultSet.ScalarFieldLists.end();
         ++resultItr) 
      resultItr->appendField(Field<Dimension, Scalar>("splat" + (*fieldItr)->name(), (*fieldItr)->nodeList()));
    for (typename vector<FieldList<Dimension, Vector> >::iterator resultItr = 
           resultSet.VectorFieldLists.begin();
         resultItr < resultSet.VectorFieldLists.end();
         ++resultItr) 
      resultItr->appendField(Field<Dimension, Vector>("splat" + (*fieldItr)->name(), (*fieldItr)->nodeList()));
    for (typename vector<FieldList<Dimension, Tensor> >::iterator resultItr = 
           resultSet.TensorFieldLists.begin();
         resultItr < resultSet.TensorFieldLists.end();
         ++resultItr) 
      resultItr->appendField(Field<Dimension, Tensor>("splat" + (*fieldItr)->name(), (*fieldItr)->nodeList()));
    for (typename vector<FieldList<Dimension, SymTensor> >::iterator resultItr = 
           resultSet.SymTensorFieldLists.begin();
         resultItr < resultSet.SymTensorFieldLists.end();
         ++resultItr) 
      resultItr->appendField(Field<Dimension, SymTensor>("splat" + (*fieldItr)->name(), (*fieldItr)->nodeList()));
  }
  vector< vector<bool> > flagNodeDone(position.numFields());
  for (typename FieldList<Dimension, Vector>::const_iterator fieldItr = position.begin();
       fieldItr < position.end(); 
       ++fieldItr) {
    flagNodeDone[fieldItr - position.begin()].resize((*fieldItr)->nodeListPtr()->numInternalNodes(), false);
  }

  // Loop over all the positions in the donor fieldList.
  for (InternalNodeIterator<Dimension> nodeItr = position.internalNodeBegin();
       nodeItr < position.internalNodeEnd();
       ++nodeItr) {

    // Check if this node has been done yet.
    if (!flagNodeDone[nodeItr.fieldID()][nodeItr.nodeID()]) {

      // Set the neighbor info over the positions we're sampling to.
      position.setMasterNodeLists(position(nodeItr), Hfield(nodeItr));
      samplePositions.setMasterNodeLists(position(nodeItr), Hfield(nodeItr));

      // Loop over the set of master nodes in the FieldList we're sampling from.
      for (MasterNodeIterator<Dimension> masterItr = position.masterNodeBegin();
           masterItr < position.masterNodeEnd();
           ++masterItr) {
        CHECK(flagNodeDone[masterItr.fieldID()][masterItr.nodeID()] == false);
   
        // Sample node (i) state.
        const Vector& ri = position(masterItr);
        const SymTensor& Hi = Hfield(masterItr);
        const Scalar& weighti = weight(masterItr);

        // Refine the set of nodes we're sampling to for this position.
        samplePositions.setRefineNodeLists(ri, Hi);

        // Loop over the refined neighbors, and determine the normalization
        // constant.
        Scalar totalWeight = 1.0e-30;
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
          totalWeight += weightij*Wij;
        }
        CHECK(totalWeight > 0.0);
        Scalar normalization = 1.0/totalWeight;

        // Loop over the refined neighbors again, and do the splat of the donor nodes
        // value to each of the sample nodes.
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

          // Loop over all the FieldLists we're sampling from, and add their contributions
          // to their correspoding result FieldList.
          const Scalar localWeight = weightij*Wij*normalization;
          for (int i = 0; i < fieldListSet.ScalarFieldLists.size(); ++i) {
            const FieldList<Dimension, Scalar>& fieldList = fieldListSet.ScalarFieldLists[i];
            FieldList<Dimension, Scalar>& result = resultSet.ScalarFieldLists[i];
            result(neighborItr) += fieldList(masterItr)*localWeight;
          }
          for (int i = 0; i < fieldListSet.VectorFieldLists.size(); ++i) {
            const FieldList<Dimension, Vector>& fieldList = fieldListSet.VectorFieldLists[i];
            FieldList<Dimension, Vector>& result = resultSet.VectorFieldLists[i];
            result(neighborItr) += fieldList(masterItr)*localWeight;
          }
          for (int i = 0; i < fieldListSet.TensorFieldLists.size(); ++i) {
            const FieldList<Dimension, Tensor>& fieldList = fieldListSet.TensorFieldLists[i];
            FieldList<Dimension, Tensor>& result = resultSet.TensorFieldLists[i];
            result(neighborItr) += fieldList(masterItr)*localWeight;
          }
          for (int i = 0; i < fieldListSet.SymTensorFieldLists.size(); ++i) {
            const FieldList<Dimension, SymTensor>& fieldList = fieldListSet.SymTensorFieldLists[i];
            FieldList<Dimension, SymTensor>& result = resultSet.SymTensorFieldLists[i];
            result(neighborItr) += fieldList(masterItr)*localWeight;
          }
        }
   
        // Flag this master node as done.
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
      cerr << "Error in FieldList::splatFieldsMash: Not all values determined on exit "
           << checkcount << endl;
    }
    CHECK(checkcount == 0);
  }

  return resultSet;
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

//============================== splatFieldsMash() ==============================
template 
FieldListSet< Dim<1> >
splatMultipleFieldsMash< Dim<1> >(const FieldListSet< Dim<1> >& fieldListSet,
                                  const FieldList<Dim<1>, Dim<1>::Vector>& position,
                                  const FieldList<Dim<1>, Dim<1>::Scalar>& weight,
                                  const FieldList<Dim<1>, Dim<1>::SymTensor>& Hfield,
                                  const TableKernel< Dim<1> >& kernel,
                                  const FieldList<Dim<1>, Dim<1>::Vector>& samplePositions,
                                  const FieldList<Dim<1>, Dim<1>::Scalar>& sampleWeight,
                                  const FieldList<Dim<1>, Dim<1>::SymTensor>& sampleHfield);

template 
FieldListSet< Dim<2> >
splatMultipleFieldsMash< Dim<2> >(const FieldListSet< Dim<2> >& fieldListSet,
                                  const FieldList<Dim<2>, Dim<2>::Vector>& position,
                                  const FieldList<Dim<2>, Dim<2>::Scalar>& weight,
                                  const FieldList<Dim<2>, Dim<2>::SymTensor>& Hfield,
                                  const TableKernel< Dim<2> >& kernel,
                                  const FieldList<Dim<2>, Dim<2>::Vector>& samplePositions,
                                  const FieldList<Dim<2>, Dim<2>::Scalar>& sampleWeight,
                                  const FieldList<Dim<2>, Dim<2>::SymTensor>& sampleHfield);

template 
FieldListSet< Dim<3> >
splatMultipleFieldsMash< Dim<3> >(const FieldListSet< Dim<3> >& fieldListSet,
                                  const FieldList<Dim<3>, Dim<3>::Vector>& position,
                                  const FieldList<Dim<3>, Dim<3>::Scalar>& weight,
                                  const FieldList<Dim<3>, Dim<3>::SymTensor>& Hfield,
                                  const TableKernel< Dim<3> >& kernel,
                                  const FieldList<Dim<3>, Dim<3>::Vector>& samplePositions,
                                  const FieldList<Dim<3>, Dim<3>::Scalar>& sampleWeight,
                                  const FieldList<Dim<3>, Dim<3>::SymTensor>& sampleHfield);

}
}
