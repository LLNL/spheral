namespace Spheral {
namespace Material {

//------------------------------------------------------------------------------
// Unit length in meters.
//------------------------------------------------------------------------------
inline
double
PhysicalConstants::unitLengthMeters() const {
  return mUnitLm;
}

//------------------------------------------------------------------------------
// Unit mass in kg.
//------------------------------------------------------------------------------
inline
double
PhysicalConstants::unitMassKg() const {
  return mUnitMkg;
}

//------------------------------------------------------------------------------
// Unit time in sec.
//------------------------------------------------------------------------------
inline
double
PhysicalConstants::unitTimeSec() const {
  return mUnitTsec;
}

//------------------------------------------------------------------------------
// Proton mass.
//------------------------------------------------------------------------------
inline
double
PhysicalConstants::protonMass() const {
  return ProtonMass;
}

//------------------------------------------------------------------------------
// Electron mass.
//------------------------------------------------------------------------------
inline
double
PhysicalConstants::electronMass() const {
  return ElectronMass;
}

//------------------------------------------------------------------------------
// Electron charge.
//------------------------------------------------------------------------------
inline
double
PhysicalConstants::electronCharge() const {
  return ElectronCharge;
}

//------------------------------------------------------------------------------
// Gravitational constant.
//------------------------------------------------------------------------------
inline
double
PhysicalConstants::G() const {
  return GGravity;
}

//------------------------------------------------------------------------------
// Speed of light.
//------------------------------------------------------------------------------
inline
double
PhysicalConstants::c() const {
  return cLight;
}

//------------------------------------------------------------------------------
// Boltzmann constant.
//------------------------------------------------------------------------------
inline
double
PhysicalConstants::kB() const {
  return kBoltzmann;
}

//------------------------------------------------------------------------------
// Avagadro's constant.
//------------------------------------------------------------------------------
inline
double
PhysicalConstants::Navogadro() const {
  return NAvogadro;
}

//------------------------------------------------------------------------------
// molar constant
//------------------------------------------------------------------------------
inline
double
PhysicalConstants::molarGasConstant() const {
  return MolarGasConstant;
}

//------------------------------------------------------------------------------
// 
//------------------------------------------------------------------------------
inline
double
PhysicalConstants::kelvinsToEnergyPerMole() const {
  return KelvinsToEnergyPerMole;
}

}
}
