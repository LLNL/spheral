#include <vector>

#include "Field/Field.hh"
#include "Geometry/Dimension.hh"

using namespace Spheral;
using namespace Spheral::FieldSpace;

typedef Spheral::FieldSpace::Field<Spheral::Dim<1>, Spheral::Dim<1>::Scalar> ScalarField1d;
typedef Spheral::FieldSpace::Field<Spheral::Dim<1>, Spheral::Dim<1>::Vector> VectorField1d;
typedef Spheral::FieldSpace::Field<Spheral::Dim<1>, Spheral::Dim<1>::Tensor> TensorField1d;
typedef Spheral::FieldSpace::Field<Spheral::Dim<1>, Spheral::Dim<1>::SymTensor> SymTensorField1d;
