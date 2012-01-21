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
// Return a monotonic smoothed estimate of the given FieldList by using a
// MASH like prescription.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
FieldList<Dimension, DataType>
smoothFieldsMash(const FieldList<Dimension, DataType>& fieldList,
                 const FieldList<Dimension, typename Dimension::Vector>& position,
                 const FieldList<Dimension, typename Dimension::Scalar>& weight,
                 const FieldList<Dimension, typename Dimension::SymTensor>& Hfield,
                 const TableKernel<Dimension>& kernel) {

  // Some convenient typedefs.
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  // Return FieldList.
  FieldList<Dimension, DataType> result;
  vector< vector<bool> > flagNodeDone(fieldList.numFields());
  result.copyFields();
  for (typename FieldList<Dimension, DataType>::const_iterator fieldItr = fieldList.begin();
       fieldItr < fieldList.end(); 
       ++fieldItr) {
    result.appendField(Field<Dimension, DataType>("smooth", (*fieldItr)->nodeList()));
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

//         // Get a pointer to the i nodes MashNodeList.
//         const MashNodeList<Dimension>* iNodeListPtr = dynamic_cast<const MashNodeList<Dimension>*>(masterItr.nodeListPtr());
//         CHECK(iNodeListPtr != 0);

        // Loop over the refined neighbors.
        const Vector& ri = position(masterItr);
        const SymTensor& Hi = Hfield(masterItr);
        const Scalar& weighti = weight(masterItr);
        const DataType& fieldi = fieldList(masterItr);

        Scalar totalWeight = 0.0;

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
          CHECK(etai >= 0.0 && etaj >= 0.0);

          // Calculate the kernel estimates for each node.
          Scalar Wi = kernel(etai, 1.0);
          Scalar Wj = kernel(etaj, 1.0);

//           // Get a pointer to the j nodes MashNodeList.
//           const MashNodeList<Dimension>* jNodeListPtr = dynamic_cast<const MashNodeList<Dimension>*>(neighborItr.nodeListPtr());
//           CHECK(jNodeListPtr != 0);

//           // Are we applying the linear correction term?
//           if (iNodeListPtr->linearCorrection()) 
//             Wi *= iNodeListPtr->a()(masterItr) +
//                   iNodeListPtr->b()(masterItr).dot(rij);
              
//           if (jNodeListPtr->linearCorrection()) 
//             Wj *= jNodeListPtr->a()(neighborItr) +
//                   jNodeListPtr->b()(neighborItr).dot(rij);
            
          // Get the symmetrized kernel weighting for this node pair.
          Scalar Wij;
          switch(neighborItr.nodeListPtr()->neighbor().neighborSearchType()) {
          case NeighborSpace::GatherScatter:
            Wij = 0.5*(Wi + Wj);
            break;

          case NeighborSpace::Gather:
            Wij = Wi;
            break;

          case NeighborSpace::Scatter:
            Wij = Wj;
            break;
          }

          // Add this nodes contribution to the master value.
          totalWeight += weightj*Wij;
          result(masterItr) += fieldj*weightj*Wij;
        }

        // This master node is finished.
        CHECK(totalWeight > 0.0);
        result(masterItr) /= totalWeight;
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
//============================== smoothFieldsMash() ==============================
namespace Spheral {
namespace FieldSpace {

using KernelSpace::TableKernel;

template 
FieldList<Dim<1>, Dim<1>::Scalar> 
smoothFieldsMash<Dim<1>, Dim<1>::Scalar>(const FieldList<Dim<1>, Dim<1>::Scalar>& fieldList,
                                     const FieldList<Dim<1>, Dim<1>::Vector>& position,
                                     const FieldList<Dim<1>, Dim<1>::Scalar>& weight,
                                     const FieldList<Dim<1>, Dim<1>::SymTensor>& Hfield,
                                     const TableKernel< Dim<1> >& kernel);
template 
FieldList<Dim<1>, Dim<1>::Vector> 
smoothFieldsMash<Dim<1>, Dim<1>::Vector>(const FieldList<Dim<1>, Dim<1>::Vector>& fieldList,
                                     const FieldList<Dim<1>, Dim<1>::Vector>& position,
                                     const FieldList<Dim<1>, Dim<1>::Scalar>& weight,
                                     const FieldList<Dim<1>, Dim<1>::SymTensor>& Hfield,
                                     const TableKernel< Dim<1> >& kernel);
template 
FieldList<Dim<1>, Dim<1>::Tensor> 
smoothFieldsMash<Dim<1>, Dim<1>::Tensor>(const FieldList<Dim<1>, Dim<1>::Tensor>& fieldList,
                                     const FieldList<Dim<1>, Dim<1>::Vector>& position,
                                     const FieldList<Dim<1>, Dim<1>::Scalar>& weight,
                                     const FieldList<Dim<1>, Dim<1>::SymTensor>& Hfield,
                                     const TableKernel< Dim<1> >& kernel);
template 
FieldList<Dim<1>, Dim<1>::SymTensor> 
smoothFieldsMash<Dim<1>, Dim<1>::SymTensor>(const FieldList<Dim<1>, Dim<1>::SymTensor>& fieldList,
                                        const FieldList<Dim<1>, Dim<1>::Vector>& position,
                                        const FieldList<Dim<1>, Dim<1>::Scalar>& weight,
                                        const FieldList<Dim<1>, Dim<1>::SymTensor>& Hfield,
                                        const TableKernel< Dim<1> >& kernel);

template 
FieldList<Dim<2>, Dim<2>::Scalar> 
smoothFieldsMash<Dim<2>, Dim<2>::Scalar>(const FieldList<Dim<2>, Dim<2>::Scalar>& fieldList,
                                     const FieldList<Dim<2>, Dim<2>::Vector>& position,
                                     const FieldList<Dim<2>, Dim<2>::Scalar>& weight,
                                     const FieldList<Dim<2>, Dim<2>::SymTensor>& Hfield,
                                     const TableKernel< Dim<2> >& kernel);
template 
FieldList<Dim<2>, Dim<2>::Vector> 
smoothFieldsMash<Dim<2>, Dim<2>::Vector>(const FieldList<Dim<2>, Dim<2>::Vector>& fieldList,
                                     const FieldList<Dim<2>, Dim<2>::Vector>& position,
                                     const FieldList<Dim<2>, Dim<2>::Scalar>& weight,
                                     const FieldList<Dim<2>, Dim<2>::SymTensor>& Hfield,
                                     const TableKernel< Dim<2> >& kernel);
template 
FieldList<Dim<2>, Dim<2>::Tensor> 
smoothFieldsMash<Dim<2>, Dim<2>::Tensor>(const FieldList<Dim<2>, Dim<2>::Tensor>& fieldList,
                                     const FieldList<Dim<2>, Dim<2>::Vector>& position,
                                     const FieldList<Dim<2>, Dim<2>::Scalar>& weight,
                                     const FieldList<Dim<2>, Dim<2>::SymTensor>& Hfield,
                                     const TableKernel< Dim<2> >& kernel);
template 
FieldList<Dim<2>, Dim<2>::SymTensor> 
smoothFieldsMash<Dim<2>, Dim<2>::SymTensor>(const FieldList<Dim<2>, Dim<2>::SymTensor>& fieldList,
                                        const FieldList<Dim<2>, Dim<2>::Vector>& position,
                                        const FieldList<Dim<2>, Dim<2>::Scalar>& weight,
                                        const FieldList<Dim<2>, Dim<2>::SymTensor>& Hfield,
                                        const TableKernel< Dim<2> >& kernel);

template 
FieldList<Dim<3>, Dim<3>::Scalar> 
smoothFieldsMash<Dim<3>, Dim<3>::Scalar>(const FieldList<Dim<3>, Dim<3>::Scalar>& fieldList,
                                     const FieldList<Dim<3>, Dim<3>::Vector>& position,
                                     const FieldList<Dim<3>, Dim<3>::Scalar>& weight,
                                     const FieldList<Dim<3>, Dim<3>::SymTensor>& Hfield,
                                     const TableKernel< Dim<3> >& kernel);
template 
FieldList<Dim<3>, Dim<3>::Vector> 
smoothFieldsMash<Dim<3>, Dim<3>::Vector>(const FieldList<Dim<3>, Dim<3>::Vector>& fieldList,
                                     const FieldList<Dim<3>, Dim<3>::Vector>& position,
                                     const FieldList<Dim<3>, Dim<3>::Scalar>& weight,
                                     const FieldList<Dim<3>, Dim<3>::SymTensor>& Hfield,
                                     const TableKernel< Dim<3> >& kernel);
template 
FieldList<Dim<3>, Dim<3>::Tensor> 
smoothFieldsMash<Dim<3>, Dim<3>::Tensor>(const FieldList<Dim<3>, Dim<3>::Tensor>& fieldList,
                                     const FieldList<Dim<3>, Dim<3>::Vector>& position,
                                     const FieldList<Dim<3>, Dim<3>::Scalar>& weight,
                                     const FieldList<Dim<3>, Dim<3>::SymTensor>& Hfield,
                                     const TableKernel< Dim<3> >& kernel);
template 
FieldList<Dim<3>, Dim<3>::SymTensor> 
smoothFieldsMash<Dim<3>, Dim<3>::SymTensor>(const FieldList<Dim<3>, Dim<3>::SymTensor>& fieldList,
                                        const FieldList<Dim<3>, Dim<3>::Vector>& position,
                                        const FieldList<Dim<3>, Dim<3>::Scalar>& weight,
                                        const FieldList<Dim<3>, Dim<3>::SymTensor>& Hfield,
                                        const TableKernel< Dim<3> >& kernel);

}
}
