#ifndef __ExtendField3d_hh__
#define __ExtendField3d_hh__

#include "Geometry/Dimension.hh"
#include "Field/NodeIterators.hh"
#include "Field/Field.hh"

#include "BPLWraps/getSetItem.hh"
#include "BPLWraps/iteratorsAsList.hh"

#include "Kernel/TableKernel.hh"

using namespace Spheral;
using namespace Spheral::FieldSpace;

typedef Spheral::Dim<3> FieldDim;

typedef Spheral::Dim<3>::Scalar Scalar;
typedef Spheral::Dim<3>::Vector Vector;
typedef Spheral::Dim<3>::Vector3d Vector3d;
typedef Spheral::Dim<3>::Tensor Tensor;
typedef Spheral::Dim<3>::SymTensor SymTensor;
typedef Spheral::Dim<3>::ThirdRankTensor TRTensor;
typedef unsigned long long ULL;

typedef Field<Spheral::Dim<3>, Scalar> ScalarField;
typedef Field<Spheral::Dim<3>, Vector> VectorField;
typedef Field<Spheral::Dim<3>, Vector3d> Vector3dField;
typedef Field<Spheral::Dim<3>, Tensor> TensorField;
typedef Field<Spheral::Dim<3>, SymTensor> SymTensorField;
typedef Field<Spheral::Dim<3>, TRTensor> ThirdRankTensorField;
typedef Field<Spheral::Dim<3>, int> IntField;
typedef Field<Spheral::Dim<3>, ULL> ULLField;
typedef Field<Spheral::Dim<3>, std::vector<double> > VectorDoubleField;

#include "ExtendField.hh"

#endif
