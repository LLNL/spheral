#include "DBC.hh"
#include "Infrastructure/SpheralFunctions.hh"

namespace Spheral {
namespace SolidMaterial {

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
SolidEquationOfState<Dimension>::
SolidEquationOfState(const double referenceDensity,
                     const double etamin,
                     const double etamax,
                     const double minimumPressure,
                     const double maximumPressure):
  Material::EquationOfState<Dimension>(minimumPressure, maximumPressure),
  mReferenceDensity(referenceDensity),
  mEtaMin(etamin),
  mEtaMax(etamax) {
  REQUIRE(distinctlyGreaterThan(mReferenceDensity, 0.0));
  REQUIRE(mEtaMin <= mEtaMax);
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
SolidEquationOfState<Dimension>::
~SolidEquationOfState() {
}

//------------------------------------------------------------------------------
// Access the reference density.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
SolidEquationOfState<Dimension>::
referenceDensity() const {
  return mReferenceDensity;
}

template<typename Dimension>
inline
void
SolidEquationOfState<Dimension>::
referenceDensity(const double x) {
  mReferenceDensity = x;
}

//------------------------------------------------------------------------------
// Access the eta min.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
SolidEquationOfState<Dimension>::
etamin() const {
  return mEtaMin;
}

template<typename Dimension>
inline
void
SolidEquationOfState<Dimension>::
etamin(const double x) {
  mEtaMin = x;
}

//------------------------------------------------------------------------------
// Access the eta max.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
SolidEquationOfState<Dimension>::
etamax() const {
  return mEtaMax;
}

template<typename Dimension>
inline
void
SolidEquationOfState<Dimension>::
etamax(const double x) {
  mEtaMax = x;
}

//------------------------------------------------------------------------------
// Compute a bounded value of eta.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
SolidEquationOfState<Dimension>::
boundedEta(const double rho) const {
  REQUIRE(distinctlyGreaterThan(mReferenceDensity, 0.0));
  return std::max(mEtaMin, std::min(mEtaMax, rho/mReferenceDensity));
}

//------------------------------------------------------------------------------
// Determine if the EOS is in a valid state.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
SolidEquationOfState<Dimension>::
valid() const {
  return (mReferenceDensity > 0.0 &&
          mEtaMin <= mEtaMax);
}

}
}
