#include "Geometry/Dimension.hh"
#include "spheralWildMagicConverters.hh"
#include "safeInv.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// 1-D bounding box.
//------------------------------------------------------------------------------
template<>
inline
Dim<1>::Box
boundingBox<Dim<1> >(const Dim<1>::Vector& xi,
                     const Dim<1>::SymTensor& Hi,
                     const Dim<1>::Scalar& kernelExtent) {
  const Dim<1>::Scalar delta = kernelExtent*safeInv(Hi.xx());
  return Dim<1>::Box(xi, delta);
}

//------------------------------------------------------------------------------
// 2-D bounding box.
//------------------------------------------------------------------------------
template<>
inline
Dim<2>::Box
boundingBox<Dim<2> >(const Dim<2>::Vector& xi,
                     const Dim<2>::SymTensor& Hi,
                     const Dim<2>::Scalar& kernelExtent) {
  typedef Dim<2>::WMVector WMVector;
  const EigenStruct<2> eigen = Hi.eigenVectors();
  return Dim<2>::Box(convertVectorToWMVector<Dim<2> >(xi),
                     convertVectorToWMVector<Dim<2> >(eigen.eigenVectors.getColumn(0)),
                     convertVectorToWMVector<Dim<2> >(eigen.eigenVectors.getColumn(1)),
                     kernelExtent / (eigen.eigenValues(0)),
                     kernelExtent / (eigen.eigenValues(1)));
}

//------------------------------------------------------------------------------
// 3-D bounding box.
//------------------------------------------------------------------------------
template<>
inline
Dim<3>::Box
boundingBox<Dim<3> >(const Dim<3>::Vector& xi,
                     const Dim<3>::SymTensor& Hi,
                     const Dim<3>::Scalar& kernelExtent) {
  typedef Dim<3>::WMVector WMVector;
  const EigenStruct<3> eigen = Hi.eigenVectors();
  return Dim<3>::Box(convertVectorToWMVector<Dim<3> >(xi),
                     convertVectorToWMVector<Dim<3> >(eigen.eigenVectors.getColumn(0)),
                     convertVectorToWMVector<Dim<3> >(eigen.eigenVectors.getColumn(1)),
                     convertVectorToWMVector<Dim<3> >(eigen.eigenVectors.getColumn(2)),
                     kernelExtent / (eigen.eigenValues(0)),
                     kernelExtent / (eigen.eigenValues(1)),
                     kernelExtent / (eigen.eigenValues(2)));
}

}
