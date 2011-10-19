#ifndef __ExtendField1d_hh__
#define __ExtendField1d_hh__

#include "Geometry/Dimension.hh"
#include "Field/NodeIterators.hh"
#include "Field/Field.hh"

#include "BPLWraps/getSetItem.hh"
#include "BPLWraps/iteratorsAsList.hh"

#include "Kernel/TableKernel.hh"

using namespace Spheral;
using namespace Spheral::FieldSpace;

typedef Spheral::Dim<1> FieldDim;

typedef Spheral::Dim<1>::Scalar Scalar;
typedef Spheral::Dim<1>::Vector Vector;
typedef Spheral::Dim<1>::Vector3d Vector3d;
typedef Spheral::Dim<1>::Tensor Tensor;
typedef Spheral::Dim<1>::SymTensor SymTensor;
typedef Spheral::Dim<1>::ThirdRankTensor TRTensor;
typedef unsigned long long ULL;

typedef Field<Spheral::Dim<1>, Scalar> ScalarField;
typedef Field<Spheral::Dim<1>, Vector> VectorField;
typedef Field<Spheral::Dim<1>, Vector3d> Vector3dField;
typedef Field<Spheral::Dim<1>, Tensor> TensorField;
typedef Field<Spheral::Dim<1>, SymTensor> SymTensorField;
typedef Field<Spheral::Dim<1>, TRTensor> ThirdRankTensorField;
typedef Field<Spheral::Dim<1>, int> IntField;
typedef Field<Spheral::Dim<1>, ULL> ULLField;
typedef Field<Spheral::Dim<1>, std::vector<double> > VectorDoubleField;

#include "ExtendField.hh"

#endif
