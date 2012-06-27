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
#include "FieldListFunctions.hh"
#include "Field/FieldList.hh"
#include "Field/Field.hh"
#include "Field/NodeIterators.hh"
#include "NodeList/NodeList.hh"
#include "Neighbor/Neighbor.hh"
#include "Kernel/TableKernel.hh"
#include "Boundary/Boundary.hh"

namespace Spheral {
namespace FieldSpace {

using namespace std;
using NodeSpace::NodeList;
using NeighborSpace::Neighbor;
using KernelSpace::TableKernel;
using BoundarySpace::Boundary;

//------------------------------------------------------------------------------
// Calculate the gradient of the divergence of a Vector FieldList.
// Explicit method that performs the div and grad operations as sequential first
// derivatives.
//------------------------------------------------------------------------------
template<typename Dimension>
FieldList<Dimension, typename Dimension::Vector>
gradDivVectorFieldList
(const FieldList<Dimension, typename Dimension::Vector>& fieldList,
 const FieldList<Dimension, typename Dimension::Vector>& position,
 const FieldList<Dimension, typename Dimension::Scalar>& weight,
 const FieldList<Dimension, typename Dimension::Scalar>& mass,
 const FieldList<Dimension, typename Dimension::Scalar>& density,
 const FieldList<Dimension, typename Dimension::SymTensor>& Hfield,
 const TableKernel<Dimension>& kernel,
 const vector<Boundary<Dimension>*>&  boundaries) {

  // Some convenient typedefs.
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  // First evaluate the divergence of the input field list.
  FieldList<Dimension, Scalar> divField = divergence(fieldList,
                                                     position,
                                                     weight,
                                                     mass,
                                                     density,
                                                     Hfield,
                                                     kernel);

  // Apply boundary conditions to the divergence.
  for (typename vector<Boundary<Dimension>*>::const_iterator bcItr = boundaries.begin();
       bcItr < boundaries.end();
       ++bcItr) {
    (*bcItr)->applyFieldListGhostBoundary(divField);
  }
  for (typename vector<Boundary<Dimension>*>::const_iterator bcItr = boundaries.begin();
       bcItr < boundaries.end();
       ++bcItr) (*bcItr)->finalizeGhostBoundary();

  // Now take the gradient of this.
  FieldList<Dimension, Vector> result = gradient(divField,
                                                 position,
                                                 weight,
                                                 mass,
                                                 density,
                                                 Hfield,
                                                 kernel);

  // That's it.
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

//============================== gradDivVectorFieldList() ==============================
template 
FieldList<Dim<1>, Dim<1>::Vector> 
gradDivVectorFieldList< Dim<1> >
(const FieldList<Dim<1>, Dim<1>::Vector>& fieldList,
 const FieldList<Dim<1>, Dim<1>::Vector>& position,
 const FieldList<Dim<1>, Dim<1>::Scalar>& weight,
 const FieldList<Dim<1>, Dim<1>::Scalar>& mass,
 const FieldList<Dim<1>, Dim<1>::Scalar>& density,
 const FieldList<Dim<1>, Dim<1>::SymTensor>& Hfield,
 const TableKernel< Dim<1> >& kernel,
 const vector<Boundary<Dim<1> >*>& boundaries);

template 
FieldList<Dim<2>, Dim<2>::Vector> 
gradDivVectorFieldList< Dim<2> >
(const FieldList<Dim<2>, Dim<2>::Vector>& fieldList,
 const FieldList<Dim<2>, Dim<2>::Vector>& position,
 const FieldList<Dim<2>, Dim<2>::Scalar>& weight,
 const FieldList<Dim<2>, Dim<2>::Scalar>& mass,
 const FieldList<Dim<2>, Dim<2>::Scalar>& density,
 const FieldList<Dim<2>, Dim<2>::SymTensor>& Hfield,
 const TableKernel< Dim<2> >& kernel,
 const vector<Boundary<Dim<2> >*>& boundaries);

template 
FieldList<Dim<3>, Dim<3>::Vector> 
gradDivVectorFieldList< Dim<3> >
(const FieldList<Dim<3>, Dim<3>::Vector>& fieldList,
 const FieldList<Dim<3>, Dim<3>::Vector>& position,
 const FieldList<Dim<3>, Dim<3>::Scalar>& weight,
 const FieldList<Dim<3>, Dim<3>::Scalar>& mass,
 const FieldList<Dim<3>, Dim<3>::Scalar>& density,
 const FieldList<Dim<3>, Dim<3>::SymTensor>& Hfield,
 const TableKernel< Dim<3> >& kernel,
 const vector<Boundary<Dim<3> >*>& boundaries);

}
}
