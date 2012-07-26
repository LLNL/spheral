//---------------------------------Spheral++----------------------------------//
// PhysicalConstants -- Choose the physical units for a given Spheral run.
//
// Created by JMO, Mon Dec  6 21:36:45 PST 1999
//----------------------------------------------------------------------------//
#include "PhysicalConstants.hh"
#include "MKSUnits.hh"
#include "CGSUnits.hh"
#include "CosmologicalUnits.hh"
#include "SolarUnits.hh"
#include "Utilities/FastMath.hh"

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
namespace Spheral {
namespace Material {

// Define the basic unit conversions in terms of the unit template type.
template<typename UnitsType> const double PhysicalConstants<UnitsType>::unitLm = UnitsType::unitLm;
template<typename UnitsType> const double PhysicalConstants<UnitsType>::unitMkg = UnitsType::unitMkg;
template<typename UnitsType> const double PhysicalConstants<UnitsType>::unitTsec = UnitsType::unitTsec;

// Define the base quantities -- we'll use the MKS SI standard.
// Source -- CRC
template<typename UnitsType> const double PhysicalConstants<UnitsType>::mpMKS = 1.672648586e-27; // kg
template<typename UnitsType> const double PhysicalConstants<UnitsType>::meMKS = 0.910953447e-30; // kg
template<typename UnitsType> const double PhysicalConstants<UnitsType>::qeMKS = 1.602176565e-19; // Coulombs
template<typename UnitsType> const double PhysicalConstants<UnitsType>::qeCGS = 4.80320425e-10;  // stat-Coulombs
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

// The electron charge is kind of funny due to the definition of the Coulomb
template<typename UnitsType> const double PhysicalConstants<UnitsType>::ElectronCharge = qeCGS / (std::sqrt(1000.0*unitMkg*FastMath::cube(100.0*unitLm))/unitTsec);

// // Explicit instantiation for each of the descendent types.
// template class PhysicalConstants<MKSUnits>;
// template class PhysicalConstants<CGSUnits>;
// template class PhysicalConstants<CosmologicalUnits>;
// template class PhysicalConstants<SolarUnits>;

}
}

