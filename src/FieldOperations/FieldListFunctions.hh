//---------------------------------Spheral++----------------------------------//
// FieldListFunctions -- A set of global functions which can be applied to
// FieldLists.
//
// Created by JMO, Wed Dec  6 21:09:29 PST 2000
//----------------------------------------------------------------------------//
#include "Geometry/MathTraits.hh"

namespace Spheral {
  namespace KernelSpace {
    template<typename Dimension> class TableKernel;
  }

  namespace FieldSpace {

    template<typename Dimension, typename DataType> class FieldList;

    // Calculate a smoothed estimate of the given FieldList.
    template<typename Dimension, typename DataType>
    FieldList<Dimension, DataType>
    smoothFields(const FieldList<Dimension, DataType>& fieldList,
                 const FieldList<Dimension, typename Dimension::Vector>& position,
                 const FieldList<Dimension, typename Dimension::Scalar>& weight,
                 const FieldList<Dimension, typename Dimension::Scalar>& mass,
                 const FieldList<Dimension, typename Dimension::Scalar>& rho,
                 const FieldList<Dimension, typename Dimension::SymTensor>& Hfield,
                 const KernelSpace::TableKernel<Dimension>& kernel);

    // Calculate the gradient of the given FieldList.
    template<typename Dimension, typename DataType>
    FieldList<Dimension, typename MathTraits<Dimension, DataType>::GradientType>
    gradient(const FieldList<Dimension, DataType>& fieldList,
             const FieldList<Dimension, typename Dimension::Vector>& position,
             const FieldList<Dimension, typename Dimension::Scalar>& weight,
             const FieldList<Dimension, typename Dimension::Scalar>& mass,
             const FieldList<Dimension, typename Dimension::Scalar>& rho,
             const FieldList<Dimension, typename Dimension::SymTensor>& Hfield,
             const KernelSpace::TableKernel<Dimension>& kernel);

    // Calculate the divergence of the given FieldList.
    template<typename Dimension, typename DataType>
    FieldList<Dimension, typename MathTraits<Dimension, DataType>::DivergenceType>
    divergence(const FieldList<Dimension, DataType>& fieldList,
               const FieldList<Dimension, typename Dimension::Vector>& position,
               const FieldList<Dimension, typename Dimension::Scalar>& weight,
               const FieldList<Dimension, typename Dimension::Scalar>& mass,
               const FieldList<Dimension, typename Dimension::Scalar>& rho,
               const FieldList<Dimension, typename Dimension::SymTensor>& Hfield,
               const KernelSpace::TableKernel<Dimension>& kernel);

    // Calculate a tensor limiter appropriate for a monotonically limited gradient.
    template<typename Dimension, typename DataType>
    FieldList<Dimension, typename Dimension::SymTensor>
    limiter(const FieldList<Dimension, DataType>& fieldList,
            const FieldList<Dimension, typename MathTraits<Dimension, DataType>::GradientType>& gradient,
            const FieldList<Dimension, typename Dimension::Vector>& position,
            const FieldList<Dimension, typename Dimension::SymTensor>& Hfield,
            const KernelSpace::TableKernel<Dimension>& kernel);

  }
}
