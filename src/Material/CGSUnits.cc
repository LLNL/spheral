//---------------------------------Spheral++----------------------------------//
// CGSUnits -- The base for the CGS unit system.
//
// Created by JMO, Fri Mar 31 17:07:41 PST 2000
//----------------------------------------------------------------------------//
#include "CGSUnits.hh"
#include "PhysicalConstants.cc"

namespace Spheral {
namespace Material {

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
CGSUnits::CGSUnits(): PhysicalConstants<CGSUnits>() {}

//------------------------------------------------------------------------------
// Set the static variables.
//------------------------------------------------------------------------------
const double CGSUnits::unitLm = 0.01; // unit length in meters
const double CGSUnits::unitMkg = 0.001; // unit mass in kg
const double CGSUnits::unitTsec = 1.0; // unit time in sec

template class PhysicalConstants<CGSUnits>;

}
}
