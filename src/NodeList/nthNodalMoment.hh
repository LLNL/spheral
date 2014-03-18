//---------------------------------Spheral++----------------------------------//
// nthNodalMoment
// nodalMoments
//
// Compute the nth moment of the local nodal distribution in \eta space:
//
//    \sum_j  (\eta_i)^n W_ij
//    -----------------------
//         \sum_j W_ij
//
// Nodal moments is specialized to simultaneously compute the non-normalized 
// zeroth and normalized first moments.
//
// Created by JMO, Mon May  9 16:20:18 PDT 2011
//----------------------------------------------------------------------------//
#ifndef __Spheral_nthNodalMoment__
#define __Spheral_nthNodalMoment__

namespace Spheral {

// Forward declarations.
namespace FieldSpace {
  template<typename Dimension, typename Value> class FieldList;
}
namespace KernelSpace {
  template<typename Dimension> class TableKernel;
}

namespace NodeSpace {

// A trait class to figure out the result type.
template<typename Dimension, unsigned n> struct MomentTraits;
template<typename Dimension> struct MomentTraits<Dimension, 0U> { typedef typename Dimension::Scalar Moment; };
template<typename Dimension> struct MomentTraits<Dimension, 1U> { typedef typename Dimension::Vector Moment; };
template<typename Dimension> struct MomentTraits<Dimension, 2U> { typedef typename Dimension::SymTensor Moment; };

// Compute a particular moment.
template<typename Dimension, typename NodeListIterator, unsigned moment>
FieldSpace::FieldList<Dimension, typename MomentTraits<Dimension, moment>::Moment>
nthNodalMoment(const NodeListIterator nodeListBegin,
               const NodeListIterator nodeListEnd,
               const KernelSpace::TableKernel<Dimension>& W,
               const bool renormalize);

// Compute the non-normalized zeroth and normalized first moment.
template<typename Dimension, typename NodeListIterator>
void
zerothAndFirstNodalMoments(const NodeListIterator nodeListBegin,
                           const NodeListIterator nodeListEnd,
                           const KernelSpace::TableKernel<Dimension>& W,
                           const bool useKernelAsGradient,
                           FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& zerothMoment,
                           FieldSpace::FieldList<Dimension, typename Dimension::Vector>& firstMoment);
}
}

#endif
