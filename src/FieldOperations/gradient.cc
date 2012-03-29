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

namespace Spheral {
namespace FieldSpace {

using namespace std;
using NodeSpace::NodeList;
using NeighborSpace::Neighbor;
using KernelSpace::TableKernel;

//------------------------------------------------------------------------------
// Calculate the gradient of a FieldList.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
FieldList<Dimension, typename MathTraits<Dimension, DataType>::GradientType>
gradient(const FieldList<Dimension, DataType>& fieldList,
         const FieldList<Dimension, typename Dimension::Vector>& position,
         const FieldList<Dimension, typename Dimension::Scalar>& weight,
         const FieldList<Dimension, typename Dimension::Scalar>& mass,
         const FieldList<Dimension, typename Dimension::Scalar>& rho,
         const FieldList<Dimension, typename Dimension::SymTensor>& Hfield,
         const TableKernel<Dimension>& kernel) {

  // Typedef's to ease typing/understandability.
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename MathTraits<Dimension, DataType>::GradientType GradientType;

  // Return FieldList.
  FieldList<Dimension, GradientType> result;
  vector< vector<bool> > flagNodeDone(fieldList.numFields());
  result.copyFields();
  for (typename FieldList<Dimension, DataType>::const_iterator
         fieldItr = fieldList.begin();
       fieldItr < fieldList.end(); 
       ++fieldItr) {
    result.appendField(Field<Dimension, GradientType>("grad", (*fieldItr)->nodeList()));
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
        const Scalar& rhoi = rho(masterItr);
        const DataType& fieldi = fieldList(masterItr);

        for (RefineNodeIterator<Dimension> neighborItr = fieldList.refineNodeBegin();
             neighborItr < fieldList.refineNodeEnd();
             ++neighborItr) {

          const Vector& rj = position(neighborItr);
          const SymTensor& Hj = Hfield(neighborItr);
          const Scalar& mj = mass(neighborItr);
          const DataType& fieldj = fieldList(neighborItr);

          const Vector rij = ri - rj;
          const Vector etai = Hi*rij;
          const Vector etaj = Hj*rij;
          const Vector etaiNorm = etai.unitVector();
          const Vector etajNorm = etaj.unitVector();

          // Get the symmetrized kernel gradient for this node pair.
          Vector gradWij;
          switch((*fieldList.begin())->nodeListPtr()->neighbor().neighborSearchType()) {
          case NeighborSpace::GatherScatter:
            gradWij = 0.5*(Hi*etaiNorm*kernel.grad(etai, Hi) + 
                           Hj*etajNorm*kernel.grad(etaj, Hj));
            break;

          case NeighborSpace::Gather:
            gradWij = Hi*etaiNorm*kernel.grad(etai, Hi);
            break;

          case NeighborSpace::Scatter:
            gradWij = Hj*etajNorm*kernel.grad(etaj, Hj);
            break;

          default:
            VERIFY2(false, "Unhandled neighbor search type.");
          }

          // Add this nodes contribution to the master value.
          result(masterItr) += mj*(fieldj - fieldi)*gradWij;
        }

        // Normalize by density.
        CHECK(rhoi > 0.0);
        result(masterItr) /= rhoi;

        // This master node is finished.
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
      cerr << "Error in FieldList::gradient: Not all values determined on exit "
           << checkcount << endl;
    }
    CHECK(checkcount == 0);
  }

  return result;
}

//------------------------------------------------------------------------------
// The limiter method.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
monotonicLimiter(const typename Dimension::Scalar& dF, 
                 const typename Dimension::Scalar& dFproj,
                 const typename Dimension::Vector& nhat,
                 const double fuzz = 1.0e-15) {
  REQUIRE(fuzz > 0.0);
  REQUIRE(fuzzyEqual(nhat.magnitude2(), 1.0));
  const double dFproj2 = dFproj*dFproj;
  return max(0.0, min(1.0, max(dF*dFproj/(dFproj2 + fuzz), fuzz/(dFproj2 + fuzz))));
}

template<typename Dimension>
inline
double
monotonicLimiter(const typename Dimension::Vector& dF,
                 const typename Dimension::Vector& dFproj,
                 const typename Dimension::Vector& nhat,
                 const double fuzz = 1.0e-15) {
  REQUIRE(fuzz > 0.0);
  REQUIRE(fuzzyEqual(nhat.magnitude2(), 1.0));
  return monotonicLimiter<Dimension>(dF.dot(nhat), dFproj.dot(nhat), nhat, fuzz);
}

//------------------------------------------------------------------------------
// Apply a monotonic limiter to a pre-computed gradient.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
FieldList<Dimension, typename Dimension::SymTensor>
limiter(const FieldList<Dimension, DataType>& fieldList,
        const FieldList<Dimension, typename MathTraits<Dimension, DataType>::GradientType>& gradient,
        const FieldList<Dimension, typename Dimension::Vector>& position,
        const FieldList<Dimension, typename Dimension::SymTensor>& Hfield,
        const TableKernel<Dimension>& kernel) {

  const int maxIterations = 10;
  const double W0 = kernel(0.0, 1.0);
  const double tiny = 1.0e-15;

  // Typedef's to ease typing/understandability.
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename MathTraits<Dimension, DataType>::GradientType GradientType;

  // Return FieldList.
  FieldList<Dimension, SymTensor> result;
  vector< vector<bool> > flagNodeDone(fieldList.numFields());
  result.copyFields();
  for (typename FieldList<Dimension, DataType>::const_iterator
         fieldItr = fieldList.begin();
       fieldItr < fieldList.end(); 
       ++fieldItr) {
    result.appendField(Field<Dimension, SymTensor>("limiter", (*fieldItr)->nodeList()));
    flagNodeDone[fieldItr - fieldList.begin()].resize((*fieldItr)->nodeListPtr()->numInternalNodes(), false);
  }

  // Loop over all the elements in the input FieldList.
  for (InternalNodeIterator<Dimension> nodeItr = fieldList.internalNodeBegin();
       nodeItr != fieldList.internalNodeEnd();
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
        const DataType& fieldi = fieldList(masterItr);
        const GradientType& gradi = gradient(masterItr);

        // Prepare the result for this node.
        SymTensor phi = SymTensor::one;

        // We iterate on the limiter.
        int iter = 0;
        Scalar phimin = 0.0;
        while (not fuzzyEqual(phimin, 1.0, 1.0e-5) and iter < maxIterations) {

          // Loop over the refined neighbors.
          Scalar weightSum = 0.0;
          SymTensor phii;
          phimin = 1.0;
          for (RefineNodeIterator<Dimension> neighborItr = fieldList.refineNodeBegin();
               neighborItr != fieldList.refineNodeEnd();
               ++neighborItr) {

            // State for node j.
            const Vector& rj = position(neighborItr);
            const SymTensor& Hj = Hfield(neighborItr);
            const DataType& fieldj = fieldList(neighborItr);

            // Compute the pair-wise limiting needed.
            const Vector rji = rj - ri;
            const DataType dF = fieldj - fieldi;
            const DataType dFproj = (phi*gradi).dot(rji);
            const Vector rhat = rji.unitVector();
            const Scalar phiij = monotonicLimiter<Dimension>(dF, dFproj, rhat);

            // Increment this iterations estimate of the limiter.
            const Vector etai = Hi*rji;
            const Scalar Wi = kernel(etai.magnitude(), 1.0)/W0;
            CHECK(Wi >= 0.0 and Wi <= 1.0);
            weightSum += Wi;
            phii += Wi*phiij/(phiij*phiij + tiny);
            phimin = min(phimin, Wi*phiij + (1.0 - Wi)*phimin);
          }

          // Increment the overall limiter.
          phi = (phii*phi).Symmetric();
        }

        // Final safety.
        phi *= phimin;
        result(masterItr) = phi;
        flagNodeDone[masterItr.fieldID()][masterItr.nodeID()] = true;
      }
    }
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

//============================== gradient() ==============================
template 
FieldList<Dim<1>, MathTraits<Dim<1>, Dim<1>::Scalar>::GradientType> 
gradient<Dim<1>, Dim<1>::Scalar>(const FieldList<Dim<1>, Dim<1>::Scalar>& fieldList,
                                     const FieldList<Dim<1>, Dim<1>::Vector>& position,
                                     const FieldList<Dim<1>, Dim<1>::Scalar>& weight,
                                     const FieldList<Dim<1>, Dim<1>::Scalar>& mass,
                                     const FieldList<Dim<1>, Dim<1>::Scalar>& density,
                                     const FieldList<Dim<1>, Dim<1>::SymTensor>& Hfield,
                                     const TableKernel< Dim<1> >& kernel);
template 
FieldList<Dim<1>, MathTraits<Dim<1>, Dim<1>::Vector>::GradientType> 
gradient<Dim<1>, Dim<1>::Vector>(const FieldList<Dim<1>, Dim<1>::Vector>& fieldList,
                                     const FieldList<Dim<1>, Dim<1>::Vector>& position,
                                     const FieldList<Dim<1>, Dim<1>::Scalar>& weight,
                                     const FieldList<Dim<1>, Dim<1>::Scalar>& mass,
                                     const FieldList<Dim<1>, Dim<1>::Scalar>& density,
                                     const FieldList<Dim<1>, Dim<1>::SymTensor>& Hfield,
                                     const TableKernel< Dim<1> >& kernel);

template 
FieldList<Dim<2>, MathTraits<Dim<2>, Dim<2>::Scalar>::GradientType> 
gradient<Dim<2>, Dim<2>::Scalar>(const FieldList<Dim<2>, Dim<2>::Scalar>& fieldList,
                                     const FieldList<Dim<2>, Dim<2>::Vector>& position,
                                     const FieldList<Dim<2>, Dim<2>::Scalar>& weight,
                                     const FieldList<Dim<2>, Dim<2>::Scalar>& mass,
                                     const FieldList<Dim<2>, Dim<2>::Scalar>& density,
                                     const FieldList<Dim<2>, Dim<2>::SymTensor>& Hfield,
                                     const TableKernel< Dim<2> >& kernel);
template 
FieldList<Dim<2>, MathTraits<Dim<2>, Dim<2>::Vector>::GradientType> 
gradient<Dim<2>, Dim<2>::Vector>(const FieldList<Dim<2>, Dim<2>::Vector>& fieldList,
                                     const FieldList<Dim<2>, Dim<2>::Vector>& position,
                                     const FieldList<Dim<2>, Dim<2>::Scalar>& weight,
                                     const FieldList<Dim<2>, Dim<2>::Scalar>& mass,
                                     const FieldList<Dim<2>, Dim<2>::Scalar>& density,
                                     const FieldList<Dim<2>, Dim<2>::SymTensor>& Hfield,
                                     const TableKernel< Dim<2> >& kernel);

template 
FieldList<Dim<3>, MathTraits<Dim<3>, Dim<3>::Scalar>::GradientType> 
gradient<Dim<3>, Dim<3>::Scalar>(const FieldList<Dim<3>, Dim<3>::Scalar>& fieldList,
                                     const FieldList<Dim<3>, Dim<3>::Vector>& position,
                                     const FieldList<Dim<3>, Dim<3>::Scalar>& weight,
                                     const FieldList<Dim<3>, Dim<3>::Scalar>& mass,
                                     const FieldList<Dim<3>, Dim<3>::Scalar>& density,
                                     const FieldList<Dim<3>, Dim<3>::SymTensor>& Hfield,
                                     const TableKernel< Dim<3> >& kernel);
template 
FieldList<Dim<3>, MathTraits<Dim<3>, Dim<3>::Vector>::GradientType> 
gradient<Dim<3>, Dim<3>::Vector>(const FieldList<Dim<3>, Dim<3>::Vector>& fieldList,
                                     const FieldList<Dim<3>, Dim<3>::Vector>& position,
                                     const FieldList<Dim<3>, Dim<3>::Scalar>& weight,
                                     const FieldList<Dim<3>, Dim<3>::Scalar>& mass,
                                     const FieldList<Dim<3>, Dim<3>::Scalar>& density,
                                     const FieldList<Dim<3>, Dim<3>::SymTensor>& Hfield,
                                     const TableKernel< Dim<3> >& kernel);

//============================== limiter() ==============================
template 
FieldList<Dim<1>, Dim<1>::SymTensor> 
limiter<Dim<1>, Dim<1>::Scalar>(const FieldList<Dim<1>, Dim<1>::Scalar>& fieldList,
                                const FieldList<Dim<1>, MathTraits<Dim<1>, Dim<1>::Scalar>::GradientType>& gradient,
                                const FieldList<Dim<1>, Dim<1>::Vector>& position,
                                const FieldList<Dim<1>, Dim<1>::SymTensor>& Hfield,
                                const TableKernel< Dim<1> >& kernel);
template 
FieldList<Dim<1>, Dim<1>::SymTensor> 
limiter<Dim<1>, Dim<1>::Vector>(const FieldList<Dim<1>, Dim<1>::Vector>& fieldList,
                                const FieldList<Dim<1>, MathTraits<Dim<1>, Dim<1>::Vector>::GradientType>& gradient,
                                const FieldList<Dim<1>, Dim<1>::Vector>& position,
                                const FieldList<Dim<1>, Dim<1>::SymTensor>& Hfield,
                                const TableKernel< Dim<1> >& kernel);
template 
FieldList<Dim<2>, Dim<2>::SymTensor> 
limiter<Dim<2>, Dim<2>::Scalar>(const FieldList<Dim<2>, Dim<2>::Scalar>& fieldList,
                                const FieldList<Dim<2>, MathTraits<Dim<2>, Dim<2>::Scalar>::GradientType>& gradient,
                                const FieldList<Dim<2>, Dim<2>::Vector>& position,
                                const FieldList<Dim<2>, Dim<2>::SymTensor>& Hfield,
                                const TableKernel< Dim<2> >& kernel);
template 
FieldList<Dim<2>, Dim<2>::SymTensor> 
limiter<Dim<2>, Dim<2>::Vector>(const FieldList<Dim<2>, Dim<2>::Vector>& fieldList,
                                const FieldList<Dim<2>, MathTraits<Dim<2>, Dim<2>::Vector>::GradientType>& gradient,
                                const FieldList<Dim<2>, Dim<2>::Vector>& position,
                                const FieldList<Dim<2>, Dim<2>::SymTensor>& Hfield,
                                const TableKernel< Dim<2> >& kernel);
template 
FieldList<Dim<3>, Dim<3>::SymTensor> 
limiter<Dim<3>, Dim<3>::Scalar>(const FieldList<Dim<3>, Dim<3>::Scalar>& fieldList,
                                const FieldList<Dim<3>, MathTraits<Dim<3>, Dim<3>::Scalar>::GradientType>& gradient,
                                const FieldList<Dim<3>, Dim<3>::Vector>& position,
                                const FieldList<Dim<3>, Dim<3>::SymTensor>& Hfield,
                                const TableKernel< Dim<3> >& kernel);
template 
FieldList<Dim<3>, Dim<3>::SymTensor> 
limiter<Dim<3>, Dim<3>::Vector>(const FieldList<Dim<3>, Dim<3>::Vector>& fieldList,
                                const FieldList<Dim<3>, MathTraits<Dim<3>, Dim<3>::Vector>::GradientType>& gradient,
                                const FieldList<Dim<3>, Dim<3>::Vector>& position,
                                const FieldList<Dim<3>, Dim<3>::SymTensor>& Hfield,
                                const TableKernel< Dim<3> >& kernel);

}
}
