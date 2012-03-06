// BPL includes.
#include "boost/python.hpp"

// Spheral++ includes.
#include "Field/FieldList.hh"
#include "Field/FieldListSet.hh"
#include "FieldOperations/FieldListFunctions.hh"
#include "FieldOperations/FieldListFunctionsMash.hh"
#include "FieldOperations/FieldListSecondDerivatives.hh"
#include "FieldOperations/PairWiseFieldListFunctions.hh"
#include "Kernel/TableKernel.hh"

#include "ExtendSampleMultipleFields2Lattice.hh"

using namespace boost::python;

using namespace Spheral;
using namespace Spheral::FieldSpace;

// Define our module.

BOOST_PYTHON_MODULE(FieldOperations) {

  // smoothFields
  def("smoothScalarFields1d", smoothFields<Dim<1>, Dim<1>::Scalar>);
  def("smoothVectorFields1d", smoothFields<Dim<1>, Dim<1>::Vector>);
  def("smoothTensorFields1d", smoothFields<Dim<1>, Dim<1>::Tensor>);
  def("smoothSymTensorFields1d", smoothFields<Dim<1>, Dim<1>::SymTensor>);

  def("smoothScalarFields2d", smoothFields<Dim<2>, Dim<2>::Scalar>);
  def("smoothVectorFields2d", smoothFields<Dim<2>, Dim<2>::Vector>);
  def("smoothTensorFields2d", smoothFields<Dim<2>, Dim<2>::Tensor>);
  def("smoothSymTensorFields2d", smoothFields<Dim<2>, Dim<2>::SymTensor>);

  def("smoothScalarFields3d", smoothFields<Dim<3>, Dim<3>::Scalar>);
  def("smoothVectorFields3d", smoothFields<Dim<3>, Dim<3>::Vector>);
  def("smoothTensorFields3d", smoothFields<Dim<3>, Dim<3>::Tensor>);
  def("smoothSymTensorFields3d", smoothFields<Dim<3>, Dim<3>::SymTensor>);

  // smoothFieldsMash
  def("smoothScalarFieldsMash1d", smoothFieldsMash<Dim<1>, Dim<1>::Scalar>);
  def("smoothVectorFieldsMash1d", smoothFieldsMash<Dim<1>, Dim<1>::Vector>);
  def("smoothTensorFieldsMash1d", smoothFieldsMash<Dim<1>, Dim<1>::Tensor>);
  def("smoothSymTensorFieldsMash1d", smoothFieldsMash<Dim<1>, Dim<1>::SymTensor>);

  def("smoothScalarFieldsMash2d", smoothFieldsMash<Dim<2>, Dim<2>::Scalar>);
  def("smoothVectorFieldsMash2d", smoothFieldsMash<Dim<2>, Dim<2>::Vector>);
  def("smoothTensorFieldsMash2d", smoothFieldsMash<Dim<2>, Dim<2>::Tensor>);
  def("smoothSymTensorFieldsMash2d", smoothFieldsMash<Dim<2>, Dim<2>::SymTensor>);

  def("smoothScalarFieldsMash3d", smoothFieldsMash<Dim<3>, Dim<3>::Scalar>);
  def("smoothVectorFieldsMash3d", smoothFieldsMash<Dim<3>, Dim<3>::Vector>);
  def("smoothTensorFieldsMash3d", smoothFieldsMash<Dim<3>, Dim<3>::Tensor>);
  def("smoothSymTensorFieldsMash3d", smoothFieldsMash<Dim<3>, Dim<3>::SymTensor>);

  // sampleFieldsMash
  def("sampleScalarFieldsMash1d", sampleFieldsMash<Dim<1>, Dim<1>::Scalar>);
  def("sampleVectorFieldsMash1d", sampleFieldsMash<Dim<1>, Dim<1>::Vector>);
  def("sampleTensorFieldsMash1d", sampleFieldsMash<Dim<1>, Dim<1>::Tensor>);
  def("sampleSymTensorFieldsMash1d", sampleFieldsMash<Dim<1>, Dim<1>::SymTensor>);

  def("sampleScalarFieldsMash2d", sampleFieldsMash<Dim<2>, Dim<2>::Scalar>);
  def("sampleVectorFieldsMash2d", sampleFieldsMash<Dim<2>, Dim<2>::Vector>);
  def("sampleTensorFieldsMash2d", sampleFieldsMash<Dim<2>, Dim<2>::Tensor>);
  def("sampleSymTensorFieldsMash2d", sampleFieldsMash<Dim<2>, Dim<2>::SymTensor>);

  def("sampleScalarFieldsMash3d", sampleFieldsMash<Dim<3>, Dim<3>::Scalar>);
  def("sampleVectorFieldsMash3d", sampleFieldsMash<Dim<3>, Dim<3>::Vector>);
  def("sampleTensorFieldsMash3d", sampleFieldsMash<Dim<3>, Dim<3>::Tensor>);
  def("sampleSymTensorFieldsMash3d", sampleFieldsMash<Dim<3>, Dim<3>::SymTensor>);

  // sampleMultipleFieldsMash
  def("sampleMultipleFieldsMash1d", sampleMultipleFieldsMash<Dim<1> >);
  def("sampleMultipleFieldsMash2d", sampleMultipleFieldsMash<Dim<2> >);
  def("sampleMultipleFieldsMash3d", sampleMultipleFieldsMash<Dim<3> >);

  // splatFieldsMash
  def("splatScalarFieldsMash1d", splatFieldsMash<Dim<1>, Dim<1>::Scalar>);
  def("splatVectorFieldsMash1d", splatFieldsMash<Dim<1>, Dim<1>::Vector>);
  def("splatTensorFieldsMash1d", splatFieldsMash<Dim<1>, Dim<1>::Tensor>);
  def("splatSymTensorFieldsMash1d", splatFieldsMash<Dim<1>, Dim<1>::SymTensor>);

  def("splatScalarFieldsMash2d", splatFieldsMash<Dim<2>, Dim<2>::Scalar>);
  def("splatVectorFieldsMash2d", splatFieldsMash<Dim<2>, Dim<2>::Vector>);
  def("splatTensorFieldsMash2d", splatFieldsMash<Dim<2>, Dim<2>::Tensor>);
  def("splatSymTensorFieldsMash2d", splatFieldsMash<Dim<2>, Dim<2>::SymTensor>);

  def("splatScalarFieldsMash3d", splatFieldsMash<Dim<3>, Dim<3>::Scalar>);
  def("splatVectorFieldsMash3d", splatFieldsMash<Dim<3>, Dim<3>::Vector>);
  def("splatTensorFieldsMash3d", splatFieldsMash<Dim<3>, Dim<3>::Tensor>);
  def("splatSymTensorFieldsMash3d", splatFieldsMash<Dim<3>, Dim<3>::SymTensor>);

  // splatMultipleFieldsMash
  def("splatMultipleFieldsMash1d", splatMultipleFieldsMash<Dim<1> >);
  def("splatMultipleFieldsMash2d", splatMultipleFieldsMash<Dim<2> >);
  def("splatMultipleFieldsMash3d", splatMultipleFieldsMash<Dim<3> >);

  // gradient
  def("gradientScalar1d", gradient<Dim<1>, Dim<1>::Scalar>);
  def("gradientVector1d", gradient<Dim<1>, Dim<1>::Vector>);

  def("gradientScalar2d", gradient<Dim<2>, Dim<2>::Scalar>);
  def("gradientVector2d", gradient<Dim<2>, Dim<2>::Vector>);

  def("gradientScalar3d", gradient<Dim<3>, Dim<3>::Scalar>);
  def("gradientVector3d", gradient<Dim<3>, Dim<3>::Vector>);

  // divergence
  def("divergenceVector1d", divergence<Dim<1>, Dim<1>::Vector>);
  def("divergenceTensor1d", divergence<Dim<1>, Dim<1>::Tensor>);
  def("divergenceSymTensor1d", divergence<Dim<1>, Dim<1>::SymTensor>);

  def("divergenceVector2d", divergence<Dim<2>, Dim<2>::Vector>);
  def("divergenceTensor2d", divergence<Dim<2>, Dim<2>::Tensor>);
  def("divergenceSymTensor2d", divergence<Dim<2>, Dim<2>::SymTensor>);

  def("divergenceVector3d", divergence<Dim<3>, Dim<3>::Vector>);
  def("divergenceTensor3d", divergence<Dim<3>, Dim<3>::Tensor>);
  def("divergenceSymTensor3d", divergence<Dim<3>, Dim<3>::SymTensor>);

  // gradientMash
  def("gradientMashScalar1d", gradientMash<Dim<1>, Dim<1>::Scalar>);
  def("gradientMashVector1d", gradientMash<Dim<1>, Dim<1>::Vector>);

  def("gradientMashScalar2d", gradientMash<Dim<2>, Dim<2>::Scalar>);
  def("gradientMashVector2d", gradientMash<Dim<2>, Dim<2>::Vector>);

  def("gradientMashScalar3d", gradientMash<Dim<3>, Dim<3>::Scalar>);
  def("gradientMashVector3d", gradientMash<Dim<3>, Dim<3>::Vector>);

  // gradientPairWise
  def("gradientPairWiseScalar1d", gradientPairWise<Dim<1>, Dim<1>::Scalar>);
  def("gradientPairWiseVector1d", gradientPairWise<Dim<1>, Dim<1>::Vector>);

  def("gradientPairWiseScalar2d", gradientPairWise<Dim<2>, Dim<2>::Scalar>);
  def("gradientPairWiseVector2d", gradientPairWise<Dim<2>, Dim<2>::Vector>);

  def("gradientPairWiseScalar3d", gradientPairWise<Dim<3>, Dim<3>::Scalar>);
  def("gradientPairWiseVector3d", gradientPairWise<Dim<3>, Dim<3>::Vector>);

  // divergencePairWise
  def("divergencePairWiseVector1d", divergencePairWise<Dim<1>, Dim<1>::Vector>);
  def("divergencePairWiseVector2d", divergencePairWise<Dim<2>, Dim<2>::Vector>);
  def("divergencePairWiseVector3d", divergencePairWise<Dim<3>, Dim<3>::Vector>);

  // Sample2Lattice (SPH)
  def("sampleMultipleFields2Lattice", sampleMultipleFields2Lattice1d);
  def("sampleMultipleFields2Lattice", sampleMultipleFields2Lattice2d);
  def("sampleMultipleFields2Lattice", sampleMultipleFields2Lattice3d);

  // Sample2Lattice (MASH)
  def("sampleMultipleFields2LatticeMash", sampleMultipleFields2LatticeMash1d);
  def("sampleMultipleFields2LatticeMash", sampleMultipleFields2LatticeMash2d);
  def("sampleMultipleFields2LatticeMash", sampleMultipleFields2LatticeMash3d);

  // limiter.
  def("limiterScalar1d", limiter<Dim<1>, Dim<1>::Scalar>);
  def("limiterVector1d", limiter<Dim<1>, Dim<1>::Vector>);
  def("limiterScalar2d", limiter<Dim<2>, Dim<2>::Scalar>);
  def("limiterVector2d", limiter<Dim<2>, Dim<2>::Vector>);
  def("limiterScalar3d", limiter<Dim<3>, Dim<3>::Scalar>);
  def("limiterVector3d", limiter<Dim<3>, Dim<3>::Vector>);
}
