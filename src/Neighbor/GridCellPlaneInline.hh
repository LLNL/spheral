//------------------------------------------------------------------------------
// Output (ostream) operator.
//------------------------------------------------------------------------------
#include "GridCellIndex.hh"
template<typename Dimension>
inline
std::ostream&
operator<<(std::ostream& os, const Spheral::NeighborSpace::GridCellPlane<Dimension>& gp) {
  os << "GridCellPlane(" << gp.point() << ", " << gp.normal() << ")";
  return os;
}
