//------------------------------------------------------------------------------
// Provide a python-centric interface to the C++ sampleMultipleFields2Lattice
// function.
//------------------------------------------------------------------------------
#ifndef __GCCXML__
#include "boost/python.hpp"
#include "boost/python/detail/api_placeholder.hpp" // This is where len lives!!!
#else
#include "../fakeboost.hh"
#endif

#include "FieldOperations/sampleMultipleFields2Lattice.hh"
#include "Geometry/Dimension.hh"

namespace Spheral {
namespace FieldSpace {

  // Forward declarations.
boost::python::tuple
sampleMultipleFields2Lattice1d(const FieldListSet<Dim<1> >& fieldListSet,
                               const FieldList<Dim<1>, Dim<1>::Vector>& position,
                               const FieldList<Dim<1>, Dim<1>::Scalar>& weight,
                               const FieldList<Dim<1>, Dim<1>::SymTensor>& Hfield,
                               const FieldList<Dim<1>, int>& mask,
                               const KernelSpace::TableKernel<Dim<1> >& W,
                               const Dim<1>::Vector& xmin,
                               const Dim<1>::Vector& xmax,
                               const boost::python::tuple& nsample0);
boost::python::tuple
sampleMultipleFields2Lattice2d(const FieldListSet<Dim<2> >& fieldListSet,
                               const FieldList<Dim<2>, Dim<2>::Vector>& position,
                               const FieldList<Dim<2>, Dim<2>::Scalar>& weight,
                               const FieldList<Dim<2>, Dim<2>::SymTensor>& Hfield,
                               const FieldList<Dim<2>, int>& mask,
                               const KernelSpace::TableKernel<Dim<2> >& W,
                               const Dim<2>::Vector& xmin,
                               const Dim<2>::Vector& xmax,
                               const boost::python::tuple& nsample0);
boost::python::tuple
sampleMultipleFields2Lattice3d(const FieldListSet<Dim<3> >& fieldListSet,
                               const FieldList<Dim<3>, Dim<3>::Vector>& position,
                               const FieldList<Dim<3>, Dim<3>::Scalar>& weight,
                               const FieldList<Dim<3>, Dim<3>::SymTensor>& Hfield,
                               const FieldList<Dim<3>, int>& mask,
                               const KernelSpace::TableKernel<Dim<3> >& W,
                               const Dim<3>::Vector& xmin,
                               const Dim<3>::Vector& xmax,
                               const boost::python::tuple& nsample0);

boost::python::tuple
sampleMultipleFields2LatticeMash1d(const FieldListSet<Dim<1> >& fieldListSet,
                                   const FieldList<Dim<1>, Dim<1>::Vector>& position,
                                   const FieldList<Dim<1>, Dim<1>::Scalar>& weight,
                                   const FieldList<Dim<1>, Dim<1>::SymTensor>& Hfield,
                                   const FieldList<Dim<1>, int>& mask,
                                   const KernelSpace::TableKernel<Dim<1> >& W,
                                   const Dim<1>::Vector& xmin,
                                   const Dim<1>::Vector& xmax,
                                   const boost::python::tuple& nsample0);
boost::python::tuple
sampleMultipleFields2LatticeMash2d(const FieldListSet<Dim<2> >& fieldListSet,
                                   const FieldList<Dim<2>, Dim<2>::Vector>& position,
                                   const FieldList<Dim<2>, Dim<2>::Scalar>& weight,
                                   const FieldList<Dim<2>, Dim<2>::SymTensor>& Hfield,
                                   const FieldList<Dim<2>, int>& mask,
                                   const KernelSpace::TableKernel<Dim<2> >& W,
                                   const Dim<2>::Vector& xmin,
                                   const Dim<2>::Vector& xmax,
                                   const boost::python::tuple& nsample0);
boost::python::tuple
sampleMultipleFields2LatticeMash3d(const FieldListSet<Dim<3> >& fieldListSet,
                                   const FieldList<Dim<3>, Dim<3>::Vector>& position,
                                   const FieldList<Dim<3>, Dim<3>::Scalar>& weight,
                                   const FieldList<Dim<3>, Dim<3>::SymTensor>& Hfield,
                                   const FieldList<Dim<3>, int>& mask,
                                   const KernelSpace::TableKernel<Dim<3> >& W,
                                   const Dim<3>::Vector& xmin,
                                   const Dim<3>::Vector& xmax,
                                   const boost::python::tuple& nsample0);

#ifndef __GCCXML__
//------------------------------------------------------------------------------
// Generic method to handle converting from boost::tuples to python::tuples.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
boost::python::tuple
wrapSampleMultipleFields(boost::tuple< std::vector< std::vector<typename Dimension::Scalar> >,
                                       std::vector< std::vector<typename Dimension::Vector> >,
                                       std::vector< std::vector<typename Dimension::Tensor> >,
                                       std::vector< std::vector<typename Dimension::SymTensor> > >
                         (*realMethod)(const FieldListSet<Dimension>& ,
                                       const FieldList<Dimension, typename Dimension::Vector>& ,
                                       const FieldList<Dimension, typename Dimension::Scalar>& ,
                                       const FieldList<Dimension, typename Dimension::SymTensor>& ,
                                       const FieldList<Dimension, int>& ,
                                       const KernelSpace::TableKernel<Dimension>& ,
                                       const typename Dimension::Vector& ,
                                       const typename Dimension::Vector& ,
                                       const std::vector<int>& ),
                         const FieldListSet<Dimension >& fieldListSet,
                         const FieldList<Dimension, typename Dimension::Vector>& position,
                         const FieldList<Dimension, typename Dimension::Scalar>& weight,
                         const FieldList<Dimension, typename Dimension::SymTensor>& Hfield,
                         const FieldList<Dimension, int>& mask,
                         const KernelSpace::TableKernel<Dimension >& W,
                         const typename Dimension::Vector& xmin,
                         const typename Dimension::Vector& xmax,
                         const boost::python::tuple& nsample0) {

  // Check the user input.
  VERIFY(boost::python::len(nsample0) == Dimension::nDim);

  // Construct an STL vector<int> from the nsample python tuple.
  std::vector<int> nsample;
  for (int i = 0; i != Dimension::nDim; ++i) nsample.push_back(boost::python::extract<int>(nsample0[i]));

  // Now call the C++ routine.
  boost::tuple< std::vector< std::vector<typename Dimension::Scalar> >,
                std::vector< std::vector<typename Dimension::Vector> >,
                std::vector< std::vector<typename Dimension::Tensor> >,
                std::vector< std::vector<typename Dimension::SymTensor> > >
  thpt = (*realMethod)(fieldListSet,
                       position,
                       weight,
                       Hfield,
                       mask,
                       W,
                       xmin,
                       xmax,
                       nsample);

  // Convert to python tuples, and include the field names.

  // Scalars.
  boost::python::list scalars;
  {
    const std::vector<std::vector<typename Dimension::Scalar> >& vecs = boost::tuples::get<0>(thpt);
    VERIFY(vecs.size() == fieldListSet.ScalarFieldLists.size());
    for (int i = 0; i != vecs.size(); ++i) {
      scalars.append(boost::python::make_tuple(fieldListSet.ScalarFieldLists[i][0]->name(),
                                               vecs[i]));
    }
  }

  // Vectors.
  boost::python::list vectors;
  {
    const std::vector<std::vector<typename Dimension::Vector> >& vecs = boost::tuples::get<1>(thpt);
    VERIFY(vecs.size() == fieldListSet.VectorFieldLists.size());
    for (int i = 0; i != vecs.size(); ++i) {
      vectors.append(boost::python::make_tuple(fieldListSet.VectorFieldLists[i][0]->name(),
                                               vecs[i]));
    }
  }

  // Tensors.
  boost::python::list tensors;
  {
    const std::vector<std::vector<typename Dimension::Tensor> >& vecs = boost::tuples::get<2>(thpt);
    VERIFY(vecs.size() == fieldListSet.TensorFieldLists.size());
    for (int i = 0; i != vecs.size(); ++i) {
      tensors.append(boost::python::make_tuple(fieldListSet.TensorFieldLists[i][0]->name(),
                                               vecs[i]));
    }
  }

  // SymTensors.
  boost::python::list symtensors;
  {
    const std::vector<std::vector<typename Dimension::SymTensor> >& vecs = boost::tuples::get<3>(thpt);
    VERIFY(vecs.size() == fieldListSet.SymTensorFieldLists.size());
    for (int i = 0; i != vecs.size(); ++i) {
      symtensors.append(boost::python::make_tuple(fieldListSet.SymTensorFieldLists[i][0]->name(),
                                                  vecs[i]));
    }
  }

  // Convert the result to a python tuple, and return it.
  VERIFY(boost::python::len(scalars) == fieldListSet.ScalarFieldLists.size());
  VERIFY(boost::python::len(vectors) == fieldListSet.VectorFieldLists.size());
  VERIFY(boost::python::len(tensors) == fieldListSet.TensorFieldLists.size());
  VERIFY(boost::python::len(symtensors) == fieldListSet.SymTensorFieldLists.size());
  return boost::python::make_tuple(scalars, vectors, tensors, symtensors);
}

//------------------------------------------------------------------------------
// Dimension specific specializations (SPH).
//------------------------------------------------------------------------------
boost::python::tuple
sampleMultipleFields2Lattice1d(const FieldListSet<Dim<1> >& fieldListSet,
                               const FieldList<Dim<1>, Dim<1>::Vector>& position,
                               const FieldList<Dim<1>, Dim<1>::Scalar>& weight,
                               const FieldList<Dim<1>, Dim<1>::SymTensor>& Hfield,
                               const FieldList<Dim<1>, int>& mask,
                               const KernelSpace::TableKernel<Dim<1> >& W,
                               const Dim<1>::Vector& xmin,
                               const Dim<1>::Vector& xmax,
                               const boost::python::tuple& nsample0) {
  return wrapSampleMultipleFields<Dim<1> >(&sampleMultipleFields2Lattice<Dim<1> >,
                                           fieldListSet,
                                           position,
                                           weight,
                                           Hfield,
                                           mask,
                                           W,
                                           xmin,
                                           xmax,
                                           nsample0);
}

boost::python::tuple
sampleMultipleFields2Lattice2d(const FieldListSet<Dim<2> >& fieldListSet,
                               const FieldList<Dim<2>, Dim<2>::Vector>& position,
                               const FieldList<Dim<2>, Dim<2>::Scalar>& weight,
                               const FieldList<Dim<2>, Dim<2>::SymTensor>& Hfield,
                               const FieldList<Dim<2>, int>& mask,
                               const KernelSpace::TableKernel<Dim<2> >& W,
                               const Dim<2>::Vector& xmin,
                               const Dim<2>::Vector& xmax,
                               const boost::python::tuple& nsample0) {
  return wrapSampleMultipleFields<Dim<2> >(&sampleMultipleFields2Lattice<Dim<2> >,
                                           fieldListSet,
                                           position,
                                           weight,
                                           Hfield,
                                           mask,
                                           W,
                                           xmin,
                                           xmax,
                                           nsample0);
}

boost::python::tuple
sampleMultipleFields2Lattice3d(const FieldListSet<Dim<3> >& fieldListSet,
                               const FieldList<Dim<3>, Dim<3>::Vector>& position,
                               const FieldList<Dim<3>, Dim<3>::Scalar>& weight,
                               const FieldList<Dim<3>, Dim<3>::SymTensor>& Hfield,
                               const FieldList<Dim<3>, int>& mask,
                               const KernelSpace::TableKernel<Dim<3> >& W,
                               const Dim<3>::Vector& xmin,
                               const Dim<3>::Vector& xmax,
                               const boost::python::tuple& nsample0) {
  return wrapSampleMultipleFields<Dim<3> >(&sampleMultipleFields2Lattice<Dim<3> >,
                                           fieldListSet,
                                           position,
                                           weight,
                                           Hfield,
                                           mask,
                                           W,
                                           xmin,
                                           xmax,
                                           nsample0);
}

//------------------------------------------------------------------------------
// Dimension specific specializations (MASH).
//------------------------------------------------------------------------------
boost::python::tuple
sampleMultipleFields2LatticeMash1d(const FieldListSet<Dim<1> >& fieldListSet,
                                   const FieldList<Dim<1>, Dim<1>::Vector>& position,
                                   const FieldList<Dim<1>, Dim<1>::Scalar>& weight,
                                   const FieldList<Dim<1>, Dim<1>::SymTensor>& Hfield,
                                   const FieldList<Dim<1>, int>& mask,
                                   const KernelSpace::TableKernel<Dim<1> >& W,
                                   const Dim<1>::Vector& xmin,
                                   const Dim<1>::Vector& xmax,
                                   const boost::python::tuple& nsample0) {
  return wrapSampleMultipleFields<Dim<1> >(&sampleMultipleFields2LatticeMash<Dim<1> >,
                                           fieldListSet,
                                           position,
                                           weight,
                                           Hfield,
                                           mask,
                                           W,
                                           xmin,
                                           xmax,
                                           nsample0);
}

boost::python::tuple
sampleMultipleFields2LatticeMash2d(const FieldListSet<Dim<2> >& fieldListSet,
                                   const FieldList<Dim<2>, Dim<2>::Vector>& position,
                                   const FieldList<Dim<2>, Dim<2>::Scalar>& weight,
                                   const FieldList<Dim<2>, Dim<2>::SymTensor>& Hfield,
                                   const FieldList<Dim<2>, int>& mask,
                                   const KernelSpace::TableKernel<Dim<2> >& W,
                                   const Dim<2>::Vector& xmin,
                                   const Dim<2>::Vector& xmax,
                                   const boost::python::tuple& nsample0) {
  return wrapSampleMultipleFields<Dim<2> >(&sampleMultipleFields2LatticeMash<Dim<2> >,
                                           fieldListSet,
                                           position,
                                           weight,
                                           Hfield,
                                           mask,
                                           W,
                                           xmin,
                                           xmax,
                                           nsample0);
}

boost::python::tuple
sampleMultipleFields2LatticeMash3d(const FieldListSet<Dim<3> >& fieldListSet,
                                   const FieldList<Dim<3>, Dim<3>::Vector>& position,
                                   const FieldList<Dim<3>, Dim<3>::Scalar>& weight,
                                   const FieldList<Dim<3>, Dim<3>::SymTensor>& Hfield,
                                   const FieldList<Dim<3>, int>& mask,
                                   const KernelSpace::TableKernel<Dim<3> >& W,
                                   const Dim<3>::Vector& xmin,
                                   const Dim<3>::Vector& xmax,
                                   const boost::python::tuple& nsample0) {
  return wrapSampleMultipleFields<Dim<3> >(&sampleMultipleFields2LatticeMash<Dim<3> >,
                                           fieldListSet,
                                           position,
                                           weight,
                                           Hfield,
                                           mask,
                                           W,
                                           xmin,
                                           xmax,
                                           nsample0);
}


#endif

}
}
