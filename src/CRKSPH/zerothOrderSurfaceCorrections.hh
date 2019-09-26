//------------------------------------------------------------------------------
// Look for any points that touch a surface (multi-material or void).
// For such points:
//   - enforce only zeroth order corrections.
//------------------------------------------------------------------------------
#include "Field/FieldList.hh"

namespace Spheral {

template<typename Dimension>
void
zerothOrderSurfaceCorrections(FieldList<Dimension, typename Dimension::Scalar>& A,
                              FieldList<Dimension, typename Dimension::Vector>& B,
                              FieldList<Dimension, typename Dimension::Tensor>& C,
                              FieldList<Dimension, typename Dimension::Vector>& gradA,
                              FieldList<Dimension, typename Dimension::Tensor>& gradB,
                              FieldList<Dimension, typename Dimension::ThirdRankTensor>& gradC,
                              const FieldList<Dimension, typename Dimension::Scalar>& m0,
                              const FieldList<Dimension, typename Dimension::Vector>& gradm0,
                              const FieldList<Dimension, int>& surfacePoint);

}
