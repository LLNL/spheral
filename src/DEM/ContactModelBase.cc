//---------------------------------Spheral++----------------------------------//
// contact Model -- root abstract class for all DEM contact models in Spheral++
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
~ContactModelBase()  {}


}
