//---------------------------------Spheral++----------------------------------//
// PhysicalConstants -- Choose the physical units for a given Spheral run.
//
// Created by JMO, Mon Dec  6 21:36:45 PST 1999
//----------------------------------------------------------------------------//
#include "PhysicalConstants.hh"
#include "Utilities/FastMath.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
PhysicalConstants::
PhysicalConstants(const double unitLm,
                  const double unitMkg,
                  const double unitTsec):
  mUnitLm(unitLm),
  mUnitMkg(unitMkg),
  mUnitTsec(unitTsec),
  ProtonMass(mpMKS/unitMkg),
  ElectronMass(meMKS/unitMkg),
  GGravity(GMKS/(unitLm/unitMkg*FastMath::square(unitLm/unitTsec))),
  cLight(cMKS/(unitLm/unitTsec)),
  kBoltzmann(kBMKS/(unitMkg*FastMath::square(unitLm/unitTsec))),
  MolarGasConstant(RgasMKS/(unitMkg*FastMath::square(unitLm/unitTsec))),
  KelvinsToEnergyPerMole(unitMkg*FastMath::square(unitLm/unitTsec)/kBMKS*NAvogadro),
  UnitMassDensity(unitMkg/unitLm/unitLm/unitLm),
  Sigma(StefanBoltzmannMKS/unitMkg*unitTsec*unitTsec*unitTsec),

  // The electron charge is kind of funny due to the definition of the Coulomb
  ElectronCharge(qeCGS / (std::sqrt(1000.0*unitMkg*FastMath::cube(100.0*unitLm))/unitTsec)) {
}

//------------------------------------------------------------------------------
// Define the base quantities -- we'll use the MKS SI standard.
// Source -- CRC
//------------------------------------------------------------------------------
const double PhysicalConstants::mpMKS =    1.672621777e-27;  // kg
const double PhysicalConstants::meMKS =     9.10938291e-31;  // kg
const double PhysicalConstants::qeMKS =     1.602176565e-19; // Coulombs
const double PhysicalConstants::qeCGS =     4.80320425e-10;  // stat-Coulombs
const double PhysicalConstants::GMKS =      6.67384e-11;     // N*m^2/kg^2
const double PhysicalConstants::cMKS =      2.99792458e8;    // m/s
const double PhysicalConstants::kBMKS =     1.3806488e-23;   // J/K
const double PhysicalConstants::RgasMKS =   8.3144621;       // J/mole/K
const double PhysicalConstants::NAvogadro = 6.02214129e23;   // mol^-1
const double PhysicalConstants::StefanBoltzmannMKS = 5.67e-8;   // W/m^2/K^4

}
