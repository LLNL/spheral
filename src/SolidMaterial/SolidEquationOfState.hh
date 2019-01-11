//---------------------------------Spheral++----------------------------------//
// SolidEquationOfState
//
// Abstract base class for the solid material equations of state.
//
// Created by JMO, Thu May  5 16:07:36 PDT 2005
//----------------------------------------------------------------------------//
#ifndef SolidEquationOfState_HH
#define SolidEquationOfState_HH

#include "Material/EquationOfState.hh"

#include <limits>

namespace Spheral {

template<typename Dimension>
class SolidEquationOfState: public EquationOfState<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  // Constructors, destructors.
  SolidEquationOfState(const double referenceDensity,
                       const double etamin,
                       const double etamax,
                       const PhysicalConstants& constants,
                       const double minimumPressure,
                       const double maximumPressure,
                       const MaterialPressureMinType minPressureType);
  virtual ~SolidEquationOfState();

  // Access the member data.
  virtual double referenceDensity() const;
  double etamin() const;
  double etamax() const;
  
  virtual void referenceDensity(const double x);
  void etamin(double x);
  void etamax(double x);
  
  // Compute eta = rho/refrho, bounded to be in [etamin, etamax].
  double boundedEta(const double rho) const;

  virtual bool valid() const;

private:
  //--------------------------- Private Interface ---------------------------//
  double mReferenceDensity, mEtaMin, mEtaMax;

  // Disallow default constructor
  SolidEquationOfState();
};

}

#include "SolidEquationOfStateInline.hh"

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class SolidEquationOfState;
}

#endif
