//---------------------------------Spheral++----------------------------------//
// MashFieldListFunctions -- A set of global functions which can be applied to
// FieldLists.
// Split out the experimental MASH prescriptions.
//
// Created by JMO, Thu Mar 28 22:15:53 PST 2002
//----------------------------------------------------------------------------//
#include <vector>
#include "Geometry/MathTraits.hh"

namespace Spheral {

template<typename Dimension> class TableKernel;
template<typename Dimension> class Boundary;

template<typename Dimension, typename DataType> class FieldList;
template<typename Dimension> class FieldListSet;

// Calculate a monotonic smoothed estimate of the given FieldList.
template<typename Dimension, typename DataType>
FieldList<Dimension, DataType>
smoothFieldsMash(const FieldList<Dimension, DataType>& fieldList,
                 const FieldList<Dimension, typename Dimension::Vector>& position,
                 const FieldList<Dimension, typename Dimension::Scalar>& weight,
                 const FieldList<Dimension, typename Dimension::SymTensor>& Hfield,
                 const TableKernel<Dimension>& kernel);

// Sample the given FieldList to a new set of positions.  Primarily useful for viz.
template<typename Dimension, typename DataType>
FieldList<Dimension, DataType>
sampleFieldsMash(const FieldList<Dimension, DataType>& fieldList,
                 const FieldList<Dimension, typename Dimension::Vector>& position,
                 const FieldList<Dimension, typename Dimension::Scalar>& weight,
                 const FieldList<Dimension, typename Dimension::SymTensor>& Hfield,
                 const TableKernel<Dimension>& kernel,
                 const FieldList<Dimension, typename Dimension::Vector>& samplePositions,
                 const FieldList<Dimension, typename Dimension::Scalar>& sampleWeight,
                 const FieldList<Dimension, typename Dimension::SymTensor>& sampleHfield);

// Same as above, but do a set of FieldLists at the same time.
template<typename Dimension>
FieldListSet<Dimension>
sampleMultipleFieldsMash(const FieldListSet<Dimension>& fieldListSet,
                         const FieldList<Dimension, typename Dimension::Vector>& position,
                         const FieldList<Dimension, typename Dimension::Scalar>& weight,
                         const FieldList<Dimension, typename Dimension::SymTensor>& Hfield,
                         const TableKernel<Dimension>& kernel,
                         const FieldList<Dimension, typename Dimension::Vector>& samplePositions,
                         const FieldList<Dimension, typename Dimension::Scalar>& sampleWeight,
                         const FieldList<Dimension, typename Dimension::SymTensor>& sampleHfield);

// Conservatively sample the given FieldList to a new set of positions.
// Primarily useful for viz.
template<typename Dimension, typename DataType>
FieldList<Dimension, DataType>
splatFieldsMash(const FieldList<Dimension, DataType>& fieldList,
                const FieldList<Dimension, typename Dimension::Vector>& position,
                const FieldList<Dimension, typename Dimension::Scalar>& weight,
                const FieldList<Dimension, typename Dimension::SymTensor>& Hfield,
                const TableKernel<Dimension>& kernel,
                const FieldList<Dimension, typename Dimension::Vector>& samplePositions,
                const FieldList<Dimension, typename Dimension::Scalar>& sampleWeight,
                const FieldList<Dimension, typename Dimension::SymTensor>& sampleHfield);

// Same as above, but do a set of FieldLists at the same time.
template<typename Dimension>
FieldListSet<Dimension>
splatMultipleFieldsMash(const FieldListSet<Dimension>& fieldListSet,
                        const FieldList<Dimension, typename Dimension::Vector>& position,
                        const FieldList<Dimension, typename Dimension::Scalar>& weight,
                        const FieldList<Dimension, typename Dimension::SymTensor>& Hfield,
                        const TableKernel<Dimension>& kernel,
                        const FieldList<Dimension, typename Dimension::Vector>& samplePositions,
                        const FieldList<Dimension, typename Dimension::Scalar>& sampleWeight,
                        const FieldList<Dimension, typename Dimension::SymTensor>& sampleHfield,
                        const std::vector<Boundary<Dimension>*>& boundaries);

// Calculate a monotonic smoothed estimate of the given FieldList.
template<typename Dimension, typename DataType>
FieldList<Dimension, DataType>
smoothFieldsMash2(const FieldList<Dimension, DataType>& fieldList,
                  const FieldList<Dimension, typename Dimension::Vector>& position,
                  const FieldList<Dimension, typename Dimension::Scalar>& weight,
                  const FieldList<Dimension, typename Dimension::Scalar>& weightDensity,
                  const FieldList<Dimension, typename Dimension::SymTensor>& Hfield,
                  const TableKernel<Dimension>& kernel);

// Calculate the gradient of the given FieldList.
template<typename Dimension, typename DataType>
FieldList<Dimension, typename MathTraits<Dimension, DataType>::GradientType>
gradientMash(const FieldList<Dimension, DataType>& fieldList,
             const FieldList<Dimension, typename Dimension::Vector>& position,
             const FieldList<Dimension, typename Dimension::Scalar>& weight,
             const FieldList<Dimension, typename Dimension::SymTensor>& Hfield,
             const TableKernel<Dimension>& kernel);

template<typename Dimension, typename DataType>
FieldList<Dimension, typename MathTraits<Dimension, DataType>::GradientType>
gradientMash2(const FieldList<Dimension, DataType>& fieldList,
              const FieldList<Dimension, typename Dimension::Vector>& position,
              const FieldList<Dimension, typename Dimension::Scalar>& weight,
              const FieldList<Dimension, typename Dimension::Scalar>& weightDensity,
              const FieldList<Dimension, typename Dimension::SymTensor>& Hfield,
              const TableKernel<Dimension>& kernel);

// Calculate the divergence of the given FieldList.
template<typename Dimension, typename DataType>
FieldList<Dimension, typename MathTraits<Dimension, DataType>::DivergenceType>
divergenceMash(const FieldList<Dimension, DataType>& fieldList,
               const FieldList<Dimension, typename Dimension::Vector>& position,
               const FieldList<Dimension, typename Dimension::Scalar>& weight,
               const FieldList<Dimension, typename Dimension::SymTensor>& Hfield,
               const TableKernel<Dimension>& kernel);

}
