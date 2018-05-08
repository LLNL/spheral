//---------------------------------Spheral++------------------------------------
// Detect surface particles leveraging the zeroth and first moments
//------------------------------------------------------------------------------
#ifndef __Spheral__detectSurface__
#define __Spheral__detectSurface__

namespace Spheral {
    // Forward declarations.
    namespace NeighborSpace {
        template<typename Dimension> class ConnectivityMap;
    }
    namespace KernelSpace {
        template<typename Dimension> class TableKernel;
    }
    namespace FieldSpace {
        template<typename Dimension, typename DataType> class FieldList;
    }
    
    namespace CRKSPHSpace {
        template<typename Dimension>
        void
        detectSurface(const NeighborSpace::ConnectivityMap<Dimension>& connectivityMap,
                      const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& m0,
                      const FieldSpace::FieldList<Dimension, typename Dimension::Vector>& m1,
                      const FieldSpace::FieldList<Dimension, typename Dimension::Vector>& position,
                      const FieldSpace::FieldList<Dimension, typename Dimension::SymTensor>& H,
                      const double detectThreshold,
                      const double detectRange,
                      const double sweepAngle,
                      FieldSpace::FieldList<Dimension, int>& surfacePoint);
    }
}

#endif
