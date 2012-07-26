//---------------------------------Spheral++----------------------------------//
// PhysicalConstants -- Choose the physical units for a given Spheral run.
//
// Created by JMO, Mon Dec  6 21:36:45 PST 1999
//----------------------------------------------------------------------------//

#ifndef __Spheral_PhysicalConstants_hh__
#define __Spheral_PhysicalConstants_hh__

namespace Spheral {
namespace Material {

template<typename UnitsType>
class PhysicalConstants {

public:
  //--------------------------- Public Interface ---------------------------//

  double protonMass() const;
  double electronMass() const;
  double electronCharge() const;
  double G() const;
  double c() const;
  double kB() const;
  double unitLengthMeters() const;
  double unitMassKg() const;
  double unitTimeSec() const;

  static const double ProtonMass;
  static const double ElectronMass;
  static const double ElectronCharge;
  static const double GGravity;
  static const double cLight;
  static const double kBoltzmann;

  static const double unitLm;
  static const double unitMkg;
  static const double unitTsec;

  static const double NAvogadro;
  static const double MolarGasConstant;
  static const double KelvinsToEnergyPerMole;

private:
  //--------------------------- Private Interface ---------------------------//
  static const double mpMKS;
  static const double meMKS;
  static const double qeMKS;
  static const double qeCGS;
  static const double GMKS;
  static const double cMKS;
  static const double kBMKS;
  static const double RgasMKS;

};

}
}

#ifndef __GCCXML__
#include "PhysicalConstantsInline.hh"
#endif

#else

// Forward declaration.
namespace Spheral {
  namespace Material {
    template<typename UnitsType> class PhysicalConstants;
  }
}

#endif

  
