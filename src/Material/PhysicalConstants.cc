//---------------------------------Spheral++----------------------------------//
// PhysicalConstants -- Choose the physical units for a given Spheral run.
//
// Created by JMO, Mon Dec  6 21:36:45 PST 1999
//----------------------------------------------------------------------------//
//#include <math.h>
#include "PhysicalConstants.hh"
#include "MKSUnits.hh"
#include "CGSUnits.hh"
#include "CosmologicalUnits.hh"
#include "Utilities/FastMath.hh"

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
namespace Spheral {
namespace Material {

template<> const double PhysicalConstants<MKSUnits>::unitLm = 1.0;                    // unit length in meters
template<> const double PhysicalConstants<MKSUnits>::unitMkg = 1.0;                   // unit mass in kg
template<> const double PhysicalConstants<MKSUnits>::unitTsec = 1.0;                  // unit time in sec

template<> const double PhysicalConstants<CGSUnits>::unitLm = 0.01;                   // unit length in meters
template<> const double PhysicalConstants<CGSUnits>::unitMkg = 0.001;                 // unit mass in kg
template<> const double PhysicalConstants<CGSUnits>::unitTsec = 1.0;                  // unit time in sec

template<> const double PhysicalConstants<CosmologicalUnits>::unitLm = 3.08567757e22; // unit length in meters
template<> const double PhysicalConstants<CosmologicalUnits>::unitMkg = 1.9891e36;    // unit mass in kg
template<> const double PhysicalConstants<CosmologicalUnits>::unitTsec = 3.155674e19; // unit time in sec

// Define the base quantities -- we'll use the MKS SI standard.
// Source -- CRC
template<typename UnitsType> const double PhysicalConstants<UnitsType>::mpMKS = 1.672648586e-27; // kg
template<typename UnitsType> const double PhysicalConstants<UnitsType>::meMKS = 0.910953447e-30; // kg
template<typename UnitsType> const double PhysicalConstants<UnitsType>::qeMKS = 1.6021775e-19;   // Coulombs
template<typename UnitsType> const double PhysicalConstants<UnitsType>::GMKS = 6.672041e-11;     // N*m^2/kg^2
template<typename UnitsType> const double PhysicalConstants<UnitsType>::cMKS = 2.99792458e8;     // m/s
template<typename UnitsType> const double PhysicalConstants<UnitsType>::kBMKS = 1.38066244e-32;  // J/K
template<typename UnitsType> const double PhysicalConstants<UnitsType>::RgasMKS = 8.31441;       // J/mole/K
template<typename UnitsType> const double PhysicalConstants<UnitsType>::NAvogadro = 6.0221367e23;

// Define all the physical constants based on their MKS values.
template<typename UnitsType> const double PhysicalConstants<UnitsType>::ProtonMass = mpMKS/unitMkg;
template<typename UnitsType> const double PhysicalConstants<UnitsType>::ElectronMass = meMKS/unitMkg;
template<typename UnitsType> const double PhysicalConstants<UnitsType>::GGravity = GMKS/(unitLm/unitMkg*FastMath::square(unitLm/unitTsec));
template<typename UnitsType> const double PhysicalConstants<UnitsType>::cLight = cMKS/(unitLm/unitTsec);
template<typename UnitsType> const double PhysicalConstants<UnitsType>::kBoltzmann = kBMKS/(unitMkg*FastMath::square(unitLm/unitTsec));
template<typename UnitsType> const double PhysicalConstants<UnitsType>::MolarGasConstant = RgasMKS/(unitMkg*FastMath::square(unitLm/unitTsec));
template<typename UnitsType> const double PhysicalConstants<UnitsType>::KelvinsToEnergyPerMole = unitMkg*FastMath::square(unitLm/unitTsec)/kBMKS*NAvogadro;

template<> const double PhysicalConstants<MKSUnits>::ProtonMass = mpMKS/unitMkg;
template<> const double PhysicalConstants<MKSUnits>::ElectronMass = meMKS/unitMkg;
template<> const double PhysicalConstants<MKSUnits>::ElectronCharge = qeMKS;
template<> const double PhysicalConstants<MKSUnits>::GGravity = GMKS/(unitLm/unitMkg*FastMath::square(unitLm/unitTsec));
template<> const double PhysicalConstants<MKSUnits>::cLight = cMKS/(unitLm/unitTsec);
template<> const double PhysicalConstants<MKSUnits>::kBoltzmann = kBMKS/(unitMkg*FastMath::square(unitLm/unitTsec));

template<> const double PhysicalConstants<CGSUnits>::ProtonMass = mpMKS/unitMkg;
template<> const double PhysicalConstants<CGSUnits>::ElectronMass = meMKS/unitMkg;
template<> const double PhysicalConstants<CGSUnits>::ElectronCharge = 10.0*cMKS*qeMKS;
template<> const double PhysicalConstants<CGSUnits>::GGravity = GMKS/(unitLm/unitMkg*FastMath::square(unitLm/unitTsec));
template<> const double PhysicalConstants<CGSUnits>::cLight = cMKS/(unitLm/unitTsec);
template<> const double PhysicalConstants<CGSUnits>::kBoltzmann = kBMKS/(unitMkg*FastMath::square(unitLm/unitTsec));

template<> const double PhysicalConstants<CosmologicalUnits>::ProtonMass = mpMKS/unitMkg;
template<> const double PhysicalConstants<CosmologicalUnits>::ElectronMass = meMKS/unitMkg;
template<> const double PhysicalConstants<CosmologicalUnits>::ElectronCharge = qeMKS;                // This is a punt for now.
template<> const double PhysicalConstants<CosmologicalUnits>::GGravity = GMKS/(unitLm/unitMkg*FastMath::square(unitLm/unitTsec));
template<> const double PhysicalConstants<CosmologicalUnits>::cLight = cMKS/(unitLm/unitTsec);
template<> const double PhysicalConstants<CosmologicalUnits>::kBoltzmann = kBMKS/(unitMkg*FastMath::square(unitLm/unitTsec));

// Explicit instantiation for each of the descendent types.
template class PhysicalConstants<MKSUnits>;
template class PhysicalConstants<CGSUnits>;
template class PhysicalConstants<CosmologicalUnits>;
}
}

