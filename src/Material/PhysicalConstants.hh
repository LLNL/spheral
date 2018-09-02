//---------------------------------Spheral++----------------------------------//
// PhysicalConstants -- Choose the physical units for a given Spheral run.
//
// Created by JMO, Mon Dec  6 21:36:45 PST 1999
//----------------------------------------------------------------------------//

#ifndef __Spheral_PhysicalConstants_hh__
#define __Spheral_PhysicalConstants_hh__

namespace Spheral {

class PhysicalConstants {

public:
  //--------------------------- Public Interface ---------------------------//
  
  // You must provide the base mass, length, and time in MKS.
  PhysicalConstants(const double unitLm,
                    const double unitMkg,
                    const double unitTsec);

  // The fundamental independent quantities.
  double unitLengthMeters() const;
  double unitMassKg() const;
  double unitTimeSec() const;

  // All the stuff we provide.
  double protonMass() const;
  double electronMass() const;
  double electronCharge() const;
  double G() const;
  double c() const;
  double kB() const;
  double Navogadro() const;
  double molarGasConstant() const;
  double kelvinsToEnergyPerMole() const;
  double unitMassDensity() const;
  double stefanBoltzmannConstant() const;

private:
  //--------------------------- Private Interface ---------------------------//
  // Independent variables.
  double mUnitLm, mUnitMkg, mUnitTsec;

  // Dependent variables.
  const double ProtonMass;
  const double ElectronMass;
  const double ElectronCharge;
  const double GGravity;
  const double cLight;
  const double kBoltzmann;
  const double MolarGasConstant;
  const double KelvinsToEnergyPerMole;
  const double UnitMassDensity;
  const double Sigma;

  // The reference MKS data we base our values on.
  static const double mpMKS;
  static const double meMKS;
  static const double qeMKS;
  static const double qeCGS;
  static const double GMKS;
  static const double cMKS;
  static const double kBMKS;
  static const double RgasMKS;
  static const double NAvogadro;
  static const double StefanBoltzmannMKS;

};

}

#include "PhysicalConstantsInline.hh"

#else

// Forward declaration.
namespace Spheral {
  class PhysicalConstants;
}

#endif

  
