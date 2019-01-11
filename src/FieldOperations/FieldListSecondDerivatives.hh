//---------------------------------Spheral++----------------------------------//
// FieldListSecondDerivatives
// A set of experimental methods for evaluating the second derivative of a 
// FieldList.
//
// Created by JMO, Wed Dec 18 22:46:54 PST 2002
//----------------------------------------------------------------------------//
#include <vector>

namespace Spheral {

template<typename Dimension> class TableKernel;
template<typename Dimension> class Boundary;

template<typename Dimension, typename DataType> class FieldList;

// Calculate the gradient of the divergence of a Vector FieldList.

// Explicit method that performs the div and grad operations as sequential first
// derivatives.
template<typename Dimension>
FieldList<Dimension, typename Dimension::Vector>
gradDivVectorFieldList
(const FieldList<Dimension, typename Dimension::Vector>& fieldList,
 const FieldList<Dimension, typename Dimension::Vector>& position,
 const FieldList<Dimension, typename Dimension::Scalar>& weight,
 const FieldList<Dimension, typename Dimension::Scalar>& mass,
 const FieldList<Dimension, typename Dimension::Scalar>& rho,
 const FieldList<Dimension, typename Dimension::SymTensor>& Hfield,
 const TableKernel<Dimension>& kernel,
 const std::vector<Boundary<Dimension>*>& boundaries);

// Explicit method that performs the div and grad operations as sequential first
// derivatives.  In this case we use the direct pairwise operators.
template<typename Dimension>
FieldList<Dimension, typename Dimension::Vector>
gradDivVectorFieldListPairWise
(const FieldList<Dimension, typename Dimension::Vector>& fieldList,
 const FieldList<Dimension, typename Dimension::Vector>& position,
 const FieldList<Dimension, typename Dimension::Scalar>& weight,
 const FieldList<Dimension, typename Dimension::Scalar>& mass,
 const FieldList<Dimension, typename Dimension::Scalar>& rho,
 const FieldList<Dimension, typename Dimension::SymTensor>& Hfield,
 const TableKernel<Dimension>& kernel);

// Simplest method, just use the second derivative of the kernel.
template<typename Dimension>
FieldList<Dimension, typename Dimension::Vector>
gradDivVectorFieldListSimple
(const FieldList<Dimension, typename Dimension::Vector>& fieldList,
 const FieldList<Dimension, typename Dimension::Vector>& position,
 const FieldList<Dimension, typename Dimension::Scalar>& weight,
 const FieldList<Dimension, typename Dimension::Scalar>& mass,
 const FieldList<Dimension, typename Dimension::Scalar>& rho,
 const FieldList<Dimension, typename Dimension::SymTensor>& Hfield,
 const TableKernel<Dimension>& kernel);

// More complex method which relies on the "golden rule", in that we move the mass
// density inside the operator.
template<typename Dimension>
FieldList<Dimension, typename Dimension::Vector>
gradDivVectorFieldListGolden
(const FieldList<Dimension, typename Dimension::Vector>& fieldList,
 const FieldList<Dimension, typename Dimension::Vector>& position,
 const FieldList<Dimension, typename Dimension::Scalar>& weight,
 const FieldList<Dimension, typename Dimension::Scalar>& mass,
 const FieldList<Dimension, typename Dimension::Scalar>& rho,
 const FieldList<Dimension, typename Dimension::SymTensor>& Hfield,
 const TableKernel<Dimension>& kernel);

// Calculate the gradient of the divergence of a Vector FieldList.
// More complex method which relies on the "golden rule", in that we move the mass
// density inside the operator.
template<typename Dimension>
FieldList<Dimension, typename Dimension::Vector>
gradDivVectorFieldListGolden2
(const FieldList<Dimension, typename Dimension::Vector>& fieldList,
 const FieldList<Dimension, typename Dimension::Vector>& position,
 const FieldList<Dimension, typename Dimension::Scalar>& weight,
 const FieldList<Dimension, typename Dimension::Scalar>& mass,
 const FieldList<Dimension, typename Dimension::Scalar>& rho,
 const FieldList<Dimension, typename Dimension::SymTensor>& Hfield,
 const TableKernel<Dimension>& kernel);

// MASH formalism for estimating grad div F.
template<typename Dimension>
FieldList<Dimension, typename Dimension::Vector>
gradDivVectorFieldListMash
(const FieldList<Dimension, typename Dimension::Vector>& fieldList,
 const FieldList<Dimension, typename Dimension::Vector>& position,
 const FieldList<Dimension, typename Dimension::Scalar>& weight,
 const FieldList<Dimension, typename Dimension::Scalar>& mass,
 const FieldList<Dimension, typename Dimension::Scalar>& rho,
 const FieldList<Dimension, typename Dimension::SymTensor>& Hfield,
 const TableKernel<Dimension>& kernel);

}
