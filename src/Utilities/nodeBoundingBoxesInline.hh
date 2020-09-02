#include "Geometry/Dimension.hh"
#include "Neighbor/Neighbor.hh"
#include "safeInv.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// The bounding box for a position and H.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
std::pair<typename Dimension::Vector, typename Dimension::Vector>
boundingBox(const typename Dimension::Vector& xi,
            const typename Dimension::SymTensor& Hi,
            const typename Dimension::Scalar& kernelExtent) {
  typedef typename Dimension::Vector Vector;
  const Vector extent = Neighbor<Dimension>::HExtent(Hi, kernelExtent);
  return std::make_pair(xi - extent, xi + extent);
}

}
