//---------------------------------Spheral++----------------------------------//
// CosmologicalUnits -- The base for the Cosmological unit system.
//
// Created by JMO, Fri Mar 31 17:07:41 PST 2000
//----------------------------------------------------------------------------//
#include "CosmologicalUnits.hh"
#include "PhysicalConstants.cc"

namespace Spheral {
namespace Material {

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
CosmologicalUnits::CosmologicalUnits(): PhysicalConstants<CosmologicalUnits>() {}

//------------------------------------------------------------------------------
// Set the static variables.
//------------------------------------------------------------------------------
const double CosmologicalUnits::unitLm = 3.08567757e22; // unit length in meters
const double CosmologicalUnits::unitMkg = 1.9891e36; // unit mass in kg
const double CosmologicalUnits::unitTsec = 3.155674e19; // unit time in sec

template class PhysicalConstants<CosmologicalUnits>;
}
}
