//---------------------------------Spheral++----------------------------------//
// Physics -- root abstract class for all DEM contact models in Spheral++
//----------------------------------------------------------------------------//
#include "ContactModelBase.hh"


namespace Spheral {

//------------------------------------------------------------------------------
// Default constructor
//------------------------------------------------------------------------------
template<typename Dimension>
ContactModelBase<Dimension>::
ContactModelBase(){}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
ContactModelBase<Dimension>::
~ContactModelBase() {}


//------------------------------------------------------------------------------
// Default No-op s
//------------------------------------------------------------------------------
template<typename Dimension>
void 
ContactModelBase<Dimension>::
force(const State<Dimension>& /*state*/,
      StateDerivatives<Dimension>& /*derivs*/) const {
}

template<typename Dimension>
void
ContactModelBase<Dimension>::
torque(const State<Dimension>& /*state*/,
       StateDerivatives<Dimension>& /*derivs*/) const {
}

}
