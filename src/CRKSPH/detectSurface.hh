//---------------------------------Spheral++------------------------------------
// Detect surface particles leveraging the zeroth and first moments
//------------------------------------------------------------------------------
#ifndef __Spheral__detectSurface__
#define __Spheral__detectSurface__

namespace Spheral {

// Forward declarations.
template<typename Dimension> class ConnectivityMap;
template<typename Dimension> class TableKernel;
template<typename Dimension, typename DataType> class FieldList;

template<typename Dimension>
void
detectSurface(const ConnectivityMap<Dimension>& connectivityMap,
              const FieldList<Dimension, typename Dimension::Scalar>& m0,
              const FieldList<Dimension, typename Dimension::Vector>& m1,
              const FieldList<Dimension, typename Dimension::Vector>& position,
              const FieldList<Dimension, typename Dimension::SymTensor>& H,
              const double detectThreshold,
              const double detectRange,
              const double sweepAngle,
              FieldList<Dimension, int>& surfacePoint);

}

#endif
