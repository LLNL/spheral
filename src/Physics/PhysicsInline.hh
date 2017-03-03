#include "DataBase/DataBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"

namespace Spheral {
namespace PhysicsSpace {

//------------------------------------------------------------------------------
// Provide iterators over over the boundary conditions.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename std::vector<BoundarySpace::Boundary<Dimension>*>::iterator
Physics<Dimension>::boundaryBegin() {
  return mBoundaryConditions.begin();
}

template<typename Dimension>
inline
typename std::vector<BoundarySpace::Boundary<Dimension>*>::iterator
Physics<Dimension>::boundaryEnd() {
  return mBoundaryConditions.end();
}

//------------------------------------------------------------------------------
// Provide const iterators over over the boundary conditions.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename std::vector<BoundarySpace::Boundary<Dimension>*>::const_iterator
Physics<Dimension>::boundaryBegin() const {
  return mBoundaryConditions.begin();
}

template<typename Dimension>
inline
typename std::vector<BoundarySpace::Boundary<Dimension>*>::const_iterator
Physics<Dimension>::boundaryEnd() const {
  return mBoundaryConditions.end();
}

//------------------------------------------------------------------------------
// Access the boundary conditions.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const std::vector<BoundarySpace::Boundary<Dimension>*>&
Physics<Dimension>::boundaryConditions() const {
  return mBoundaryConditions;
}

}
}
