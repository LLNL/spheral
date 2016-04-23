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
#include "NodeList/FluidNodeList.hh"
#include "Neighbor/Neighbor.hh"
#include "Kernel/TableKernel.hh"
#include "Geometry/MathTraits.hh"

#include "Utilities/DBC.hh"

#ifdef USE_MPI
#include "Neighbor/NestedGridNeighbor.hh"
#include "Distributed/NestedGridDistributedBoundary.hh"
#include "DataBase/DataBase.hh"
#endif

namespace Spheral {
namespace FieldSpace {

using namespace std;
using NodeSpace::NodeList;
using NodeSpace::FluidNodeList;
using NeighborSpace::Neighbor;
using KernelSpace::TableKernel;
#ifdef USE_MPI
using NeighborSpace::NestedGridNeighbor;
using BoundarySpace::NestedGridDistributedBoundary;
using DataBaseSpace::DataBase;
#endif

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

#ifdef USE_MPI
  // Check what sort of DistributedBoundary we should use.
  VERIFY2(dynamic_cast<const NestedGridNeighbor<Dimension>*>(&(position[0]->nodeListPtr()->neighbor())) != 0,
          "splatMultipleFieldsMash ERROR: currently only understand how to do parallel with NestedGridNeighbor objects.");
  NestedGridDistributedBoundary<Dimension>& distributedBoundary = NestedGridDistributedBoundary<Dimension>::instance();

  // Build distributed ghost nodes.
  DataBase<Dimension> db;
  for (typename FieldList<Dimension, Vector>::const_iterator fieldItr = position.begin();
       fieldItr < position.end(); 
       ++fieldItr) {
    FluidNodeList<Dimension>& nodeList = dynamic_cast<FluidNodeList<Dimension>&>(const_cast<NodeList<Dimension>&>((*fieldItr)->nodeList()));
    db.appendNodeList(nodeList);
    nodeList.numGhostNodes(0);
    nodeList.neighbor().updateNodes();
  }
  for (typename FieldList<Dimension, Vector>::const_iterator fieldItr = samplePositions.begin();
       fieldItr < samplePositions.end(); 
       ++fieldItr) {
    FluidNodeList<Dimension>& nodeList = dynamic_cast<FluidNodeList<Dimension>&>(const_cast<NodeList<Dimension>&>((*fieldItr)->nodeList()));
    db.appendNodeList(nodeList);
    nodeList.numGhostNodes(0);
    nodeList.neighbor().updateNodes();
  }
  distributedBoundary.setAllGhostNodes(db);
  distributedBoundary.applyFieldListGhostBoundary(const_cast<FieldList<Dimension, Scalar>&>(sampleWeight));
  distributedBoundary.finalizeGhostBoundary();
  for (typename FieldList<Dimension, Vector>::const_iterator fieldItr = position.begin();
       fieldItr < position.end(); 
       ++fieldItr) {
    const_cast<NodeList<Dimension>&>((*fieldItr)->nodeList()).neighbor().updateNodes();
  }
  for (typename FieldList<Dimension, Vector>::const_iterator fieldItr = samplePositions.begin();
       fieldItr < samplePositions.end(); 
       ++fieldItr) {
    const_cast<NodeList<Dimension>&>((*fieldItr)->nodeList()).neighbor().updateNodes();
  }

  // BLAGO
  cerr << "Master nodes:" << endl;
  for (unsigned i = position[0]->nodeList().firstGhostNode(); i != position[0]->nodeList().numNodes(); ++i) cerr << "   " << position(0,i) << endl;
  cerr << " Slave nodes:" << endl;
  for (GhostNodeIterator<Dimension> itr = samplePositions.ghostNodeBegin(); itr != samplePositions.ghostNodeEnd(); ++itr) cerr << "   " << samplePositions(itr) << endl;
  // BLAGO
#endif

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
  
  FieldList<Dimension, int> flagNodeDone(FieldSpace::Copy);
  FieldList<Dimension, Scalar> normalization(FieldSpace::Copy);
  for (typename FieldList<Dimension, Vector>::const_iterator fieldItr = position.begin();
       fieldItr < position.end(); 
       ++fieldItr) {
    flagNodeDone.appendNewField("flag nodes", (*fieldItr)->nodeList(), 0);
    normalization.appendNewField("normalization", (*fieldItr)->nodeList(), 1.0e-30);
  }

  // This first pass counts the total weight from each donor node, so we can compute the normalization correctly.
  // Loop over all the positions in the donor fieldList.
  for (InternalNodeIterator<Dimension> nodeItr = position.internalNodeBegin();
       nodeItr < position.internalNodeEnd();
       ++nodeItr) {

    // Check if this node has been done yet.
    if (flagNodeDone(nodeItr) == 0) {

      // Set the neighbor info over the positions we're sampling to.
      position.setMasterNodeLists(position(nodeItr), Hfield(nodeItr));
      samplePositions.setMasterNodeLists(position(nodeItr), Hfield(nodeItr));

      // Loop over the set of master nodes in the FieldList we're sampling from.
      for (MasterNodeIterator<Dimension> masterItr = position.masterNodeBegin();
           masterItr < position.masterNodeEnd();
           ++masterItr) {
        CHECK(flagNodeDone(masterItr) == 0);
   
        // Sample node (i) state.
        const Vector& ri = position(masterItr);
        const SymTensor& Hi = Hfield(masterItr);
        const Scalar& weighti = weight(masterItr);

        // Refine the set of nodes we're sampling to for this position.
        samplePositions.setRefineNodeLists(ri, Hi);

        // Loop over the refined neighbors, and determine the normalization
        // constant.
        // cerr << position(masterItr) << " : ";
        unsigned blago = 0;
        for (RefineNodeIterator<Dimension> neighborItr = samplePositions.refineNodeBegin();
             neighborItr < samplePositions.refineNodeEnd();
             ++neighborItr) {
          // cerr << samplePositions(neighborItr) << " ";
          ++blago;

          // Node j's state.
          const Vector& rj = samplePositions(neighborItr);
          const SymTensor& Hj = sampleHfield(neighborItr);
          const Scalar& weightj = sampleWeight(neighborItr);

          const Vector rij = ri - rj;
          const Scalar etai = (Hi*rij).magnitude();
          const Scalar etaj = (Hj*rij).magnitude();
          CHECK2(etai >= 0.0 && etaj >= 0.0, etai << " " << etaj);

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
          normalization(masterItr) += weightij*Wij;
        }
        // cerr << blago << endl;
        CHECK(normalization(masterItr) > 0.0);
        normalization(masterItr) = 1.0/normalization(masterItr);
        // cerr << " NORM : " << position(masterItr) << " " << normalization(masterItr) << endl;

        // Flag this master node as done.
        flagNodeDone(masterItr) = 1;
      }
    }
  }

  // After we're done, all nodes in all NodeLists should be flagged as done.
  CHECK(flagNodeDone.min() == 1);

#ifdef USE_MPI
  // Apply boundaries to the donor information.
  distributedBoundary.applyFieldListGhostBoundary(normalization);
  for (int i = 0; i < fieldListSet.ScalarFieldLists.size(); ++i) {
    distributedBoundary.applyFieldListGhostBoundary(const_cast<FieldList<Dimension, Scalar>&>(fieldListSet.ScalarFieldLists[i]));
  }
  for (int i = 0; i < fieldListSet.VectorFieldLists.size(); ++i) {
    distributedBoundary.applyFieldListGhostBoundary(const_cast<FieldList<Dimension, Vector>&>(fieldListSet.VectorFieldLists[i]));
  }
  for (int i = 0; i < fieldListSet.TensorFieldLists.size(); ++i) {
    distributedBoundary.applyFieldListGhostBoundary(const_cast<FieldList<Dimension, Tensor>&>(fieldListSet.TensorFieldLists[i]));
  }
  for (int i = 0; i < fieldListSet.SymTensorFieldLists.size(); ++i) {
    distributedBoundary.applyFieldListGhostBoundary(const_cast<FieldList<Dimension, SymTensor>&>(fieldListSet.SymTensorFieldLists[i]));
  }
  distributedBoundary.finalizeGhostBoundary();
#endif

  FieldList<Dimension, int> flagSampleDone(FieldSpace::Copy);
  for (typename FieldList<Dimension, Vector>::const_iterator fieldItr = samplePositions.begin();
       fieldItr < samplePositions.end(); 
       ++fieldItr) {
    flagSampleDone.appendNewField("flag sample nodes", (*fieldItr)->nodeList(), 0);
  }

  // Now do a pass over the sampling nodes and splat from the donors.
  for (InternalNodeIterator<Dimension> nodeItr = samplePositions.internalNodeBegin();
       nodeItr < samplePositions.internalNodeEnd();
       ++nodeItr) {

    // Check if this node has been done yet.
    if (flagSampleDone(nodeItr) == 0) {

      // Set the neighbor info over the positions we're sampling to.
      position.setMasterNodeLists(samplePositions(nodeItr), sampleHfield(nodeItr));
      samplePositions.setMasterNodeLists(samplePositions(nodeItr), sampleHfield(nodeItr));

      // Loop over the set of master nodes in the FieldList we're sampling from.
      for (MasterNodeIterator<Dimension> masterItr = samplePositions.masterNodeBegin();
           masterItr < samplePositions.masterNodeEnd();
           ++masterItr) {
        CHECK(flagSampleDone(masterItr) == 0);
   
        // Sample node (i) state.
        const Vector& ri = samplePositions(masterItr);
        const SymTensor& Hi = sampleHfield(masterItr);
        const Scalar& weighti = sampleWeight(masterItr);

        // Refine the set of nodes we're donating from to for this position.
        position.setRefineNodeLists(ri, Hi);

        // Loop over the refined neighbors again, and do the splat of the donor nodes
        // value to each of the sample nodes.
        for (RefineNodeIterator<Dimension> neighborItr = position.refineNodeBegin();
             neighborItr < position.refineNodeEnd();
             ++neighborItr) {

          // Node j's state.
          const Vector& rj = position(neighborItr);
          const SymTensor& Hj = Hfield(neighborItr);
          const Scalar& weightj = weight(neighborItr);

          const Vector rij = ri - rj;
          const Scalar etai = (Hi*rij).magnitude();
          const Scalar etaj = (Hj*rij).magnitude();
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
            Wij = Wj;
            weightij = weightj;
            break;
         
          case NeighborSpace::Scatter:
            Wij = Wi;
            weightij = weighti;
            break;

          default:
            VERIFY2(false, "Unhandled neighbor search type.");
          }

          // Loop over all the FieldLists we're sampling from, and add their contributions
          // to their correspoding result FieldList.
          const Scalar localWeight = weightij*Wij*normalization(neighborItr);
          for (int i = 0; i < fieldListSet.ScalarFieldLists.size(); ++i) {
            const FieldList<Dimension, Scalar>& fieldList = fieldListSet.ScalarFieldLists[i];
            FieldList<Dimension, Scalar>& result = resultSet.ScalarFieldLists[i];
            result(masterItr) += fieldList(neighborItr)*localWeight;
          }
          for (int i = 0; i < fieldListSet.VectorFieldLists.size(); ++i) {
            const FieldList<Dimension, Vector>& fieldList = fieldListSet.VectorFieldLists[i];
            FieldList<Dimension, Vector>& result = resultSet.VectorFieldLists[i];
            result(masterItr) += fieldList(neighborItr)*localWeight;
          }
          for (int i = 0; i < fieldListSet.TensorFieldLists.size(); ++i) {
            const FieldList<Dimension, Tensor>& fieldList = fieldListSet.TensorFieldLists[i];
            FieldList<Dimension, Tensor>& result = resultSet.TensorFieldLists[i];
            result(masterItr) += fieldList(neighborItr)*localWeight;
          }
          for (int i = 0; i < fieldListSet.SymTensorFieldLists.size(); ++i) {
            const FieldList<Dimension, SymTensor>& fieldList = fieldListSet.SymTensorFieldLists[i];
            FieldList<Dimension, SymTensor>& result = resultSet.SymTensorFieldLists[i];
            result(masterItr) += fieldList(neighborItr)*localWeight;
          }
        }
   
        // Flag this master node as done.
        flagSampleDone(masterItr) = 1;
      }
    }
  }

  // After we're done, all nodes in all NodeLists should be flagged as done.
  CHECK(flagSampleDone.min() == 1);

  return resultSet;
}

}
}

