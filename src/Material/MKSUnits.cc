//---------------------------------Spheral++----------------------------------//
// MKSUnits -- The base for the MKS unit system.
//
// Created by JMO, Fri Mar 31 17:07:41 PST 2000
//----------------------------------------------------------------------------//
#include "MKSUnits.hh"
#include "PhysicalConstants.cc"

namespace Spheral {
namespace Material {

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
MKSUnits::MKSUnits(): PhysicalConstants<MKSUnits>() {}

//------------------------------------------------------------------------------
// Set the static variables.
//------------------------------------------------------------------------------
const double MKSUnits::unitLm = 1.0; // unit length in meters
const double MKSUnits::unitMkg = 1.0; // unit mass in kg
const double MKSUnits::unitTsec = 1.0; // unit time in sec

template<>                   const double PhysicalConstants<MKSUnits>::ElectronCharge = qeMKS;
template class PhysicalConstants<MKSUnits>;
}
}

