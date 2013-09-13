//---------------------------------Spheral++----------------------------------//
// SolidEquationOfState
//
// Abstract base class for the solid material equations of state.
//
// Created by JMO, Thu May  5 16:07:36 PDT 2005
//----------------------------------------------------------------------------//
#ifndef SolidEquationOfState_HH
#define SolidEquationOfState_HH

#include <limits>
#include "Material/EquationOfState.hh"

namespace Spheral {
namespace SolidMaterial {

template<typename Dimension>
class SolidEquationOfState: public Material::EquationOfState<Dimension> {

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
                       const Material::PhysicalConstants& constants,
                       const double minimumPressure = -std::numeric_limits<double>::max(),
                       const double maximumPressure = std::numeric_limits<double>::max());
  virtual ~SolidEquationOfState();

  // Access the member data.
  virtual double referenceDensity() const;
  double etamin() const;
  double etamax() const;
  
  virtual void referenceDensity(const double x);
  void etamin(const double x);
  void etamax(const double x);
  
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
}

#ifndef __GCCXML__
#include "SolidEquationOfStateInline.hh"
#endif

#else

// Forward declaration.
namespace Spheral {
  namespace SolidMaterial {
    template<typename Dimension> class SolidEquationOfState;
  }
}

#endif
