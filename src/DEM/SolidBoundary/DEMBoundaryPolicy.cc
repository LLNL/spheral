//---------------------------------Spheral++----------------------------------//
// DEMBoundaryPolicy -- general policy to allow solid boundaries to move in
//                      DEM. 
//
// J.M. Pearl 2023
//----------------------------------------------------------------------------//

#include "DEM/SolidBoundary/DEMBoundaryPolicy.hh"
#include "DEM/SolidBoundary/SolidBoundary.hh"


namespace Spheral {

//------------------------------------------------------------------------------
// Constructors.
//------------------------------------------------------------------------------
template<typename Dimension>
DEMBoundaryPolicy<Dimension>::
DEMBoundaryPolicy():
  UpdatePolicyBase<Dimension>(),
  mSolidBoundaries() {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
DEMBoundaryPolicy<Dimension>::
~DEMBoundaryPolicy() {
}

//------------------------------------------------------------------------------
// Update the field.
template<typename Dimension>
void
DEMBoundaryPolicy<Dimension>::
update(const KeyType& /*key*/,
       State<Dimension>& /*state*/,
       StateDerivatives<Dimension>& /*derivs*/,
       const double multiplier,
       const double t,
       const double dt) {

  for (ConstSolidBoundaryIterator solidboundItr = mSolidBoundaries.begin();
        solidboundItr != mSolidBoundaries.end();
        ++solidboundItr){
    (*solidboundItr)->update(multiplier,t,dt);
  }
}


//------------------------------------------------------------------------------
// Add a Boundary condition to the end of the current boundary list.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DEMBoundaryPolicy<Dimension>::
appendSolidBoundary(SolidBoundary<Dimension>& boundary) {
    mSolidBoundaries.push_back(&boundary);
}

//------------------------------------------------------------------------------
// Add a Boundary condition to the beginning of the current boundary list.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DEMBoundaryPolicy<Dimension>::
prependSolidBoundary(SolidBoundary<Dimension>& boundary) {
    mSolidBoundaries.insert(mSolidBoundaries.begin(), &boundary);
}

//------------------------------------------------------------------------------
// Clear (erase) the boundary condition list.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DEMBoundaryPolicy<Dimension>::
clearSolidBoundaries() {
  mSolidBoundaries = std::vector<SolidBoundary<Dimension>*>();
}

//------------------------------------------------------------------------------
// Test if the given Boundary condition is listed in the physics package.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
DEMBoundaryPolicy<Dimension>::
haveSolidBoundary(const SolidBoundary<Dimension>& boundary) const {
  return std::count(mSolidBoundaries.begin(), mSolidBoundaries.end(), &boundary) > 0;
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
DEMBoundaryPolicy<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const DEMBoundaryPolicy<Dimension>* rhsPtr = dynamic_cast<const DEMBoundaryPolicy<Dimension>*>(&rhs);
  return rhsPtr != 0;
}

}

