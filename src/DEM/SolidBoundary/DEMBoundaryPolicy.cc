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
DEMBoundaryPolicy(const std::vector<SolidBoundary<Dimension>*>& solidBoundaries):
  UpdatePolicyBase<Dimension>(),
  mSolidBoundariesRef(solidBoundaries) {
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
//------------------------------------------------------------------------------
template<typename Dimension>
void
DEMBoundaryPolicy<Dimension>::
update(const KeyType& /*key*/,
       State<Dimension>& /*state*/,
       StateDerivatives<Dimension>& /*derivs*/,
       const double multiplier,
       const double t,
       const double dt) {

  for (ConstSolidBoundaryIterator solidboundItr = mSolidBoundariesRef.begin();
        solidboundItr != mSolidBoundariesRef.end();
        ++solidboundItr){
    (*solidboundItr)->update(multiplier,t,dt);
  }
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

