//---------------------------------Spheral++----------------------------------//
// sampleMultipleFields2Lattice.
//
// SPH sample all the Fields in a FieldListSet to a lattice of positions.
// The results are returned as a 
// tuple< vector< vector<Scalar> >,
//        vector< vector<Vector> >,
//        vector< vector<Tensor> >,
//        vector< vector<SymTensor> >
//
// Created by JMO, Wed Nov 16 10:40:07 PST 2005
//----------------------------------------------------------------------------//
#include <vector>
#include "boost/tuple/tuple.hpp"

namespace Spheral {

  namespace KernelSpace {
    template<typename Dimension> class TableKernel;
  }

  namespace FieldSpace {

    template<typename Dimension, typename DataType> class FieldList;
    template<typename Dimension> class FieldListSet;

    // Simultaneously SPH sample multiple FieldLists to a lattice.
    template<typename Dimension>
    void
    sampleMultipleFields2Lattice(const FieldListSet<Dimension>& fieldListSet,
                                 const FieldList<Dimension, typename Dimension::Vector>& position,
                                 const FieldList<Dimension, typename Dimension::Scalar>& weight,
                                 const FieldList<Dimension, typename Dimension::SymTensor>& Hfield,
                                 const FieldList<Dimension, int>& mask,
                                 const KernelSpace::TableKernel<Dimension>& W,
                                 const typename Dimension::Vector& xmin,
                                 const typename Dimension::Vector& xmax,
                                 const std::vector<int>& nsample,
                                 std::vector< std::vector<typename Dimension::Scalar> >& scalarValues,
                                 std::vector< std::vector<typename Dimension::Vector> >& vectorValues,
                                 std::vector< std::vector<typename Dimension::Tensor> >& tensorValues,
                                 std::vector< std::vector<typename Dimension::SymTensor> >& symTensorValues);

    // Simultaneously MASH sample multiple FieldLists to a lattice.
    template<typename Dimension>
    void
    sampleMultipleFields2LatticeMash(const FieldListSet<Dimension>& fieldListSet,
                                     const FieldList<Dimension, typename Dimension::Vector>& position,
                                     const FieldList<Dimension, typename Dimension::Scalar>& weight,
                                     const FieldList<Dimension, typename Dimension::SymTensor>& Hfield,
                                     const FieldList<Dimension, int>& mask,
                                     const KernelSpace::TableKernel<Dimension>& W,
                                     const typename Dimension::Vector& xmin,
                                     const typename Dimension::Vector& xmax,
                                     const std::vector<int>& nsample,
                                     std::vector< std::vector<typename Dimension::Scalar> >& scalarValues,
                                     std::vector< std::vector<typename Dimension::Vector> >& vectorValues,
                                     std::vector< std::vector<typename Dimension::Tensor> >& tensorValues,
                                     std::vector< std::vector<typename Dimension::SymTensor> >& symTensorValues);

  }
}
