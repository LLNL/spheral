#ifndef __ExtendField2d_hh__
#define __ExtendField2d_hh__

#include "Geometry/Dimension.hh"
#include "Field/NodeIterators.hh"
#include "Field/Field.hh"

#include "BPLWraps/getSetItem.hh"
#include "BPLWraps/iteratorsAsList.hh"

#include "Kernel/TableKernel.hh"

using namespace Spheral;
using namespace Spheral::FieldSpace;

typedef Spheral::Dim<2> FieldDim;

typedef Spheral::Dim<2>::Scalar Scalar;
typedef Spheral::Dim<2>::Vector Vector;
typedef Spheral::Dim<2>::Vector3d Vector3d;
typedef Spheral::Dim<2>::Tensor Tensor;
typedef Spheral::Dim<2>::SymTensor SymTensor;
typedef Spheral::Dim<2>::ThirdRankTensor TRTensor;
typedef unsigned long long ULL;

typedef Field<Spheral::Dim<2>, Scalar> ScalarField;
typedef Field<Spheral::Dim<2>, Vector> VectorField;
typedef Field<Spheral::Dim<2>, Vector3d> Vector3dField;
typedef Field<Spheral::Dim<2>, Tensor> TensorField;
typedef Field<Spheral::Dim<2>, SymTensor> SymTensorField;
typedef Field<Spheral::Dim<2>, TRTensor> ThirdRankTensorField;
typedef Field<Spheral::Dim<2>, int> IntField;
typedef Field<Spheral::Dim<2>, ULL> ULLField;
typedef Field<Spheral::Dim<2>, std::vector<double> > VectorDoubleField;

#include "ExtendField.hh"

#endif
