//---------------------------------Spheral++----------------------------------//
// SolarUnits -- The base for the Solar unit system.
//
// Created by JMO, Tue Mar 20 10:35:56 PDT 2012
//----------------------------------------------------------------------------//

#include "SolarUnits.hh"

//------------------------------------------------------------------------------
// Set the static variables.
//------------------------------------------------------------------------------
namespace Spheral {
  namespace Material {
    const double SolarUnits::unitLm = 149597870700.0;   // unit length in meters (1 AU)
    const double SolarUnits::unitMkg = 1.98892e30;      // unit mass in kg (1 Msun)
    const double SolarUnits::unitTsec = 365.25*3600*24; // unit time in sec (1 year)
  }
}
