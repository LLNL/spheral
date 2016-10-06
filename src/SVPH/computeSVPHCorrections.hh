//---------------------------------Spheral++------------------------------------
// Compute the SVPH corrections.
// Based on the reproducing kernel ideas, similar to CRKSPH.
//------------------------------------------------------------------------------
#ifndef __Spheral__computeSVPHCorrections__
#define __Spheral__computeSVPHCorrections__

#include "Geometry/Dimension.hh"

namespace Spheral {

  // Forward declarations.
  namespace NeighborSpace {
    template<typename Dimension> class ConnectivityMap;
  }
  namespace KernelSpace {
    template<typename Dimension> class TableKernel;
  }
  namespace FieldSpace {
    template<typename Dimension, typename DataType> class Field;
    template<typename Dimension, typename DataType> class FieldList;
  }

  namespace SVPHSpace {

    // Single NodeList version.
    template<typename Dimension>
    void
    computeSVPHCorrections(const NeighborSpace::ConnectivityMap<Dimension>& connectivityMap,
                           const KernelSpace::TableKernel<Dimension>& W,
                           const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& volume,
                           const FieldSpace::FieldList<Dimension, typename Dimension::Vector>& position,
                           const FieldSpace::FieldList<Dimension, typename Dimension::SymTensor>& H,
                           FieldSpace::Field<Dimension, typename Dimension::Scalar>& A,
                           FieldSpace::Field<Dimension, typename Dimension::Vector>& B,
                           FieldSpace::Field<Dimension, typename Dimension::Tensor>& gradB);

    // Full FieldList version.
    template<typename Dimension>
    void
    computeSVPHCorrections(const NeighborSpace::ConnectivityMap<Dimension>& connectivityMap,
                           const KernelSpace::TableKernel<Dimension>& W,
                           const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& volume,
                           const FieldSpace::FieldList<Dimension, typename Dimension::Vector>& position,
                           const FieldSpace::FieldList<Dimension, typename Dimension::SymTensor>& H,
                           FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& A,
                           FieldSpace::FieldList<Dimension, typename Dimension::Vector>& B,
                           FieldSpace::FieldList<Dimension, typename Dimension::Tensor>& gradB);

    // // Specializations.
    // template<>
    // void
    // computeSVPHCorrections<Dim<2> >(const NeighborSpace::ConnectivityMap<Dim<2> >& connectivityMap,
    //                                 const KernelSpace::TableKernel<Dim<2> >& W,
    //                                 const FieldSpace::FieldList<Dim<2> ,  Dim<2>::Scalar>& volume,
    //                                 const FieldSpace::FieldList<Dim<2> ,  Dim<2>::Vector>& position,
    //                                 const FieldSpace::FieldList<Dim<2> ,  Dim<2>::SymTensor>& H,
    //                                 FieldSpace::FieldList<Dim<2> ,  Dim<2>::Scalar>& A,
    //                                 FieldSpace::FieldList<Dim<2> ,  Dim<2>::Vector>& B,
    //                                 FieldSpace::FieldList<Dim<2> ,  Dim<2>::Tensor>& gradB);
  }
}

#endif
