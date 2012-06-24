namespace Spheral {
namespace Material {

//------------------------------------------------------------------------------
// Proton mass.
//------------------------------------------------------------------------------
template<typename UnitsType>
inline
double
PhysicalConstants<UnitsType>::protonMass() const {
  return ProtonMass;
}

//------------------------------------------------------------------------------
// Electron mass.
//------------------------------------------------------------------------------
template<typename UnitsType>
inline
double
PhysicalConstants<UnitsType>::electronMass() const {
  return ElectronMass;
}

//------------------------------------------------------------------------------
// Electron charge.
//------------------------------------------------------------------------------
template<typename UnitsType>
inline
double
PhysicalConstants<UnitsType>::electronCharge() const {
  return ElectronCharge;
}

//------------------------------------------------------------------------------
// Gravitational constant.
//------------------------------------------------------------------------------
template<typename UnitsType>
inline
double
PhysicalConstants<UnitsType>::G() const {
  return GGravity;
}

//------------------------------------------------------------------------------
// Speed of light.
//------------------------------------------------------------------------------
template<typename UnitsType>
inline
double
PhysicalConstants<UnitsType>::c() const {
  return cLight;
}

//------------------------------------------------------------------------------
// Speed of light.
//------------------------------------------------------------------------------
template<typename UnitsType>
inline
double
PhysicalConstants<UnitsType>::kB() const {
  return kBoltzmann;
}

//------------------------------------------------------------------------------
// Unit length in meters.
//------------------------------------------------------------------------------
template<typename UnitsType>
inline
double
PhysicalConstants<UnitsType>::unitLengthMeters() const {
  return unitLm;
}

//------------------------------------------------------------------------------
// Unit mass in kg.
//------------------------------------------------------------------------------
template<typename UnitsType>
inline
double
PhysicalConstants<UnitsType>::unitMassKg() const {
  return unitMkg;
}

//------------------------------------------------------------------------------
// Unit time in sec.
//------------------------------------------------------------------------------
template<typename UnitsType>
inline
double
PhysicalConstants<UnitsType>::unitTimeSec() const {
  return unitTsec;
}

}
}
