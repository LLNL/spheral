//---------------------------------Spheral++----------------------------------//
// SteinbergGuinanLundStrength -- Implements the Steinberg, Guinan, & Lund
// rate dependent strength model.
//
//   Steinberg, D.J. & Lund, C.M. (1988), J. App. Phys., 65, 1528.
//
// Created by JMO, Wed Sep 8 15:18:44 2004
//----------------------------------------------------------------------------//
#ifndef __Spheral_SteinbergGuinanLundStrength_hh__
#define __Spheral_SteinbergGuinanLundStrength_hh__

#include "StrengthModel.hh"
#include "PolynomialFit.hh"
#include "SolidEquationOfState.hh"
#include "SteinbergGuinanStrength.hh"

namespace Spheral {

template<typename Dimension>
class SteinbergGuinanLundStrength: public SteinbergGuinanStrength<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  // Constructors, destructor.
  SteinbergGuinanLundStrength(const SolidEquationOfState<Dimension>& eos,
                              const double G0,
                              const double A,
                              const double B,
                              const double Y0,
                              const double Ymax,
                              const double Yp,
                              const double beta,
                              const double gamma0,
                              const double nhard,
                              const double C1,
                              const double C2,
                              const double UK,
                              const double YP,
                              const double YTmax,
                              const NinthOrderPolynomialFit& coldEnergyFit,
                              const NinthOrderPolynomialFit& meltEnergyFit);
  virtual ~SteinbergGuinanLundStrength();

  // Override the required generic interface.
  virtual double yieldStrength(const double density,
                               const double specificThermalEnergy,
                               const double pressure,
                               const double plasticStrain,
                               const double plasticStrainRate) const;

  // Access the strength parameters.
  double C1() const;
  double C2() const;
  double UK() const;
  double YP() const;
  double YTmax() const;

private:
  //--------------------------- Private Interface ---------------------------//
  double mC1;
  double mC2;
  double mUK;
  double mYP;
  double mYTmax;

  // No copying or assignment.
  SteinbergGuinanLundStrength(const SteinbergGuinanLundStrength&);
  SteinbergGuinanLundStrength& operator=(const SteinbergGuinanLundStrength&);
};

}

#endif

