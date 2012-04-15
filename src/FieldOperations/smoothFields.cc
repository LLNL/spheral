//---------------------------------Spheral++----------------------------------//
// FieldListFunctions -- A set of global functions which can be applied to
// FieldLists.
//
// Created by JMO, Wed Dec  6 21:09:29 PST 2000
//----------------------------------------------------------------------------//
#include "FieldListFunctions.hh"
#include "Field/FieldList.hh"
#include "Field/Field.hh"
#include "Field/NodeIterators.hh"
#include "NodeList/NodeList.hh"
#include "Neighbor/Neighbor.hh"
#include "Kernel/TableKernel.hh"
#include "Geometry/MathTraits.hh"

//------------------------------------------------------------------------------
// Return a smoothed estimate of the given FieldList.
//------------------------------------------------------------------------------
namespace Spheral {
namespace FieldSpace {

using namespace std;
using NodeSpace::NodeList;
using NeighborSpace::Neighbor;
using KernelSpace::TableKernel;

template<typename Dimension, typename DataType>
FieldList<Dimension, DataType>
smoothFields(const FieldList<Dimension, DataType>& fieldList,
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
  FieldList<Dimension, DataType> result;
  vector< vector<bool> > flagNodeDone(fieldList.numFields());
  result.copyFields();
  for (typename FieldList<Dimension, DataType>::const_iterator fieldItr = fieldList.begin();
       fieldItr < fieldList.end(); 
       ++fieldItr) {
    result.appendField(Field<Dimension, DataType>("smooth " + (*fieldItr)->name(), (*fieldItr)->nodeList()));
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

          // Get the symmetrized kernel weighting for this node pair.
          Scalar Wij;
          switch((*fieldList.begin())->nodeListPtr()->neighbor().neighborSearchType()) {
          case NeighborSpace::GatherScatter:
            Wij = 0.5*(kernel(etai, Hi) + kernel(etaj, Hj));
            break;

          case NeighborSpace::Gather:
            Wij = kernel(etai, Hi);
            break;

          case NeighborSpace::Scatter:
            Wij = kernel(etaj, Hj);
            break;

          default:
            VERIFY2(false, "Unhandled neighbor search type.");
          }

          // Add this nodes contribution to the master value.
          result(masterItr) += weightj*fieldj*Wij;
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

}
}

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
//============================== smoothFields() ==============================
namespace Spheral {
namespace FieldSpace {

using KernelSpace::TableKernel;

template 
FieldList<Dim<1>, Dim<1>::Scalar> 
smoothFields<Dim<1>, Dim<1>::Scalar>(const FieldList<Dim<1>, Dim<1>::Scalar>& fieldList,
                                     const FieldList<Dim<1>, Dim<1>::Vector>& position,
                                     const FieldList<Dim<1>, Dim<1>::Scalar>& weight,
                                     const FieldList<Dim<1>, Dim<1>::Scalar>& mass,
                                     const FieldList<Dim<1>, Dim<1>::Scalar>& density,
                                     const FieldList<Dim<1>, Dim<1>::SymTensor>& Hfield,
                                     const TableKernel< Dim<1> >& kernel);
template 
FieldList<Dim<1>, Dim<1>::Vector> 
smoothFields<Dim<1>, Dim<1>::Vector>(const FieldList<Dim<1>, Dim<1>::Vector>& fieldList,
                                     const FieldList<Dim<1>, Dim<1>::Vector>& position,
                                     const FieldList<Dim<1>, Dim<1>::Scalar>& weight,
                                     const FieldList<Dim<1>, Dim<1>::Scalar>& mass,
                                     const FieldList<Dim<1>, Dim<1>::Scalar>& density,
                                     const FieldList<Dim<1>, Dim<1>::SymTensor>& Hfield,
                                     const TableKernel< Dim<1> >& kernel);
template 
FieldList<Dim<1>, Dim<1>::Tensor> 
smoothFields<Dim<1>, Dim<1>::Tensor>(const FieldList<Dim<1>, Dim<1>::Tensor>& fieldList,
                                     const FieldList<Dim<1>, Dim<1>::Vector>& position,
                                     const FieldList<Dim<1>, Dim<1>::Scalar>& weight,
                                     const FieldList<Dim<1>, Dim<1>::Scalar>& mass,
                                     const FieldList<Dim<1>, Dim<1>::Scalar>& density,
                                     const FieldList<Dim<1>, Dim<1>::SymTensor>& Hfield,
                                     const TableKernel< Dim<1> >& kernel);
template 
FieldList<Dim<1>, Dim<1>::SymTensor> 
smoothFields<Dim<1>, Dim<1>::SymTensor>(const FieldList<Dim<1>, Dim<1>::SymTensor>& fieldList,
                                        const FieldList<Dim<1>, Dim<1>::Vector>& position,
                                        const FieldList<Dim<1>, Dim<1>::Scalar>& weight,
                                        const FieldList<Dim<1>, Dim<1>::Scalar>& mass,
                                        const FieldList<Dim<1>, Dim<1>::Scalar>& density,
                                        const FieldList<Dim<1>, Dim<1>::SymTensor>& Hfield,
                                        const TableKernel< Dim<1> >& kernel);

template 
FieldList<Dim<2>, Dim<2>::Scalar> 
smoothFields<Dim<2>, Dim<2>::Scalar>(const FieldList<Dim<2>, Dim<2>::Scalar>& fieldList,
                                     const FieldList<Dim<2>, Dim<2>::Vector>& position,
                                     const FieldList<Dim<2>, Dim<2>::Scalar>& weight,
                                     const FieldList<Dim<2>, Dim<2>::Scalar>& mass,
                                     const FieldList<Dim<2>, Dim<2>::Scalar>& density,
                                     const FieldList<Dim<2>, Dim<2>::SymTensor>& Hfield,
                                     const TableKernel< Dim<2> >& kernel);
template 
FieldList<Dim<2>, Dim<2>::Vector> 
smoothFields<Dim<2>, Dim<2>::Vector>(const FieldList<Dim<2>, Dim<2>::Vector>& fieldList,
                                     const FieldList<Dim<2>, Dim<2>::Vector>& position,
                                     const FieldList<Dim<2>, Dim<2>::Scalar>& weight,
                                     const FieldList<Dim<2>, Dim<2>::Scalar>& mass,
                                     const FieldList<Dim<2>, Dim<2>::Scalar>& density,
                                     const FieldList<Dim<2>, Dim<2>::SymTensor>& Hfield,
                                     const TableKernel< Dim<2> >& kernel);
template 
FieldList<Dim<2>, Dim<2>::Tensor> 
smoothFields<Dim<2>, Dim<2>::Tensor>(const FieldList<Dim<2>, Dim<2>::Tensor>& fieldList,
                                     const FieldList<Dim<2>, Dim<2>::Vector>& position,
                                     const FieldList<Dim<2>, Dim<2>::Scalar>& weight,
                                     const FieldList<Dim<2>, Dim<2>::Scalar>& mass,
                                     const FieldList<Dim<2>, Dim<2>::Scalar>& density,
                                     const FieldList<Dim<2>, Dim<2>::SymTensor>& Hfield,
                                     const TableKernel< Dim<2> >& kernel);
template 
FieldList<Dim<2>, Dim<2>::SymTensor> 
smoothFields<Dim<2>, Dim<2>::SymTensor>(const FieldList<Dim<2>, Dim<2>::SymTensor>& fieldList,
                                        const FieldList<Dim<2>, Dim<2>::Vector>& position,
                                        const FieldList<Dim<2>, Dim<2>::Scalar>& weight,
                                        const FieldList<Dim<2>, Dim<2>::Scalar>& mass,
                                        const FieldList<Dim<2>, Dim<2>::Scalar>& density,
                                        const FieldList<Dim<2>, Dim<2>::SymTensor>& Hfield,
                                        const TableKernel< Dim<2> >& kernel);

template 
FieldList<Dim<3>, Dim<3>::Scalar> 
smoothFields<Dim<3>, Dim<3>::Scalar>(const FieldList<Dim<3>, Dim<3>::Scalar>& fieldList,
                                     const FieldList<Dim<3>, Dim<3>::Vector>& position,
                                     const FieldList<Dim<3>, Dim<3>::Scalar>& weight,
                                     const FieldList<Dim<3>, Dim<3>::Scalar>& mass,
                                     const FieldList<Dim<3>, Dim<3>::Scalar>& density,
                                     const FieldList<Dim<3>, Dim<3>::SymTensor>& Hfield,
                                     const TableKernel< Dim<3> >& kernel);
template 
FieldList<Dim<3>, Dim<3>::Vector> 
smoothFields<Dim<3>, Dim<3>::Vector>(const FieldList<Dim<3>, Dim<3>::Vector>& fieldList,
                                     const FieldList<Dim<3>, Dim<3>::Vector>& position,
                                     const FieldList<Dim<3>, Dim<3>::Scalar>& weight,
                                     const FieldList<Dim<3>, Dim<3>::Scalar>& mass,
                                     const FieldList<Dim<3>, Dim<3>::Scalar>& density,
                                     const FieldList<Dim<3>, Dim<3>::SymTensor>& Hfield,
                                     const TableKernel< Dim<3> >& kernel);
template 
FieldList<Dim<3>, Dim<3>::Tensor> 
smoothFields<Dim<3>, Dim<3>::Tensor>(const FieldList<Dim<3>, Dim<3>::Tensor>& fieldList,
                                     const FieldList<Dim<3>, Dim<3>::Vector>& position,
                                     const FieldList<Dim<3>, Dim<3>::Scalar>& weight,
                                     const FieldList<Dim<3>, Dim<3>::Scalar>& mass,
                                     const FieldList<Dim<3>, Dim<3>::Scalar>& density,
                                     const FieldList<Dim<3>, Dim<3>::SymTensor>& Hfield,
                                     const TableKernel< Dim<3> >& kernel);
template 
FieldList<Dim<3>, Dim<3>::SymTensor> 
smoothFields<Dim<3>, Dim<3>::SymTensor>(const FieldList<Dim<3>, Dim<3>::SymTensor>& fieldList,
                                        const FieldList<Dim<3>, Dim<3>::Vector>& position,
                                        const FieldList<Dim<3>, Dim<3>::Scalar>& weight,
                                        const FieldList<Dim<3>, Dim<3>::Scalar>& mass,
                                        const FieldList<Dim<3>, Dim<3>::Scalar>& density,
                                        const FieldList<Dim<3>, Dim<3>::SymTensor>& Hfield,
                                        const TableKernel< Dim<3> >& kernel);
}
}
