//---------------------------------Spheral++----------------------------------//
// SteinbergGuinanStrength -- Implements the Steinberg-Guinan strength model.
//   ** Need reference **
//
// Created by JMO, Wed Sep 8 15:18:44 2004
//----------------------------------------------------------------------------//
#ifndef __Spheral_SteinbergGuinanStrength_hh__
#define __Spheral_SteinbergGuinanStrength_hh__

#include "StrengthModel.hh"
#ifndef __GCCXML__
#include "PolynomialFit.hh"
#endif

namespace Spheral {
  namespace SolidMaterial {
    class NinthOrderPolynomialFit;
    template<typename Dimension> class SolidEquationOfState;
  }
}

namespace Spheral {
namespace SolidMaterial {

template<typename Dimension>
class SteinbergGuinanStrength: public StrengthModel<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  // Constructors, destructor.
  SteinbergGuinanStrength(const SolidEquationOfState<Dimension>& eos,
                          const double G0,
                          const double A,
                          const double B,
                          const double Y0,
                          const double Ymax,
                          const double Yp,
                          const double beta,
                          const double gamma0,
                          const double nhard,
                          const NinthOrderPolynomialFit& coldEnergyFit,
                          const NinthOrderPolynomialFit& meltEnergyFit);
  virtual ~SteinbergGuinanStrength();

  // Override the required generic interface.
  virtual double shearModulus(const double density,
                              const double specificThermalEnergy,
                              const double pressure) const;

  virtual double yieldStrength(const double density,
                               const double specificThermalEnergy,
                               const double pressure,
                               const double plasticStrain,
                               const double plasticStrainRate) const;

  // Override the method to compute the full sound speed.
  virtual double soundSpeed(const double density,
                            const double specificThermalEnergy,
                            const double pressure,
                            const double fluidSoundSpeed) const;

  // Melt attenuation.
  double meltAttenuation(const double rho, const double eps) const;

  // Steinberg-Guinan "temperature".
  double computeTemperature(const double rho, const double eps) const;

  // Access the strength parameters.
  double G0() const;
  double A() const;
  double B() const;
  double Y0() const;
  double Ymax() const;
  double Yp() const;
  double beta() const;
  double gamma0() const;
  double nhard() const;
  const NinthOrderPolynomialFit& coldEnergyFit() const;
  const NinthOrderPolynomialFit& meltEnergyFit() const;

private:
  //--------------------------- Private Interface ---------------------------//
#ifndef __GCCXML__
  const SolidEquationOfState<Dimension>* mEOSPtr;
  double mG0;
  double mA;
  double mB;
  double mY0;
  double mYmax;
  double mYp;
  double mbeta;
  double mgamma0;
  double mnhard;
  NinthOrderPolynomialFit mColdEnergyFit;
  NinthOrderPolynomialFit mMeltEnergyFit;
#endif

  // No copying or assignment.
  SteinbergGuinanStrength(const SteinbergGuinanStrength&);
  SteinbergGuinanStrength& operator=(const SteinbergGuinanStrength&);
};

}
}

#else

// Forward declaration.
namespace Spheral {
  namespace SolidMaterial {
    template<typename Dimension> class SteinbergGuinanStrength;
  }
}

#endif

