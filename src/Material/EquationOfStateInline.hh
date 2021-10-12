namespace Spheral {

//------------------------------------------------------------------------------
// Optional molecular weight.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
EquationOfState<Dimension>::molecularWeight() const {
  VERIFY2(false, "EquationOfState " << this << " does not provide molecularWeight.");
  return 0.0;
}

//------------------------------------------------------------------------------
// Units.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const PhysicalConstants&
EquationOfState<Dimension>::constants() const {
  return mConstants;
}

//------------------------------------------------------------------------------
// Min pressure.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
EquationOfState<Dimension>::minimumPressure() const {
  return mMinimumPressure;
}

template<typename Dimension>
inline
void
EquationOfState<Dimension>::minimumPressure(const double val) {
  mMinimumPressure = val;
}

//------------------------------------------------------------------------------
// Max pressure.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
EquationOfState<Dimension>::maximumPressure() const {
  return mMaximumPressure;
}

template<typename Dimension>
inline
void
EquationOfState<Dimension>::maximumPressure(const double val) {
  mMaximumPressure = val;
}

//------------------------------------------------------------------------------
// Min damage pressure.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
EquationOfState<Dimension>::minimumDamagePressure() const {
  return mMinimumDamagePressure;
}

template<typename Dimension>
inline
void
EquationOfState<Dimension>::minimumDamagePressure(const double val) {
  mMinimumDamagePressure = val;
}

//------------------------------------------------------------------------------
// Min pressure algorithm choice.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
MaterialPressureMinType
EquationOfState<Dimension>::minimumPressureType() const {
  return mMinPressureType;
}

template<typename Dimension>
inline
void
EquationOfState<Dimension>::minimumPressureType(const MaterialPressureMinType x) {
  mMinPressureType = x;
}

//------------------------------------------------------------------------------
// Apply pressure limits to a candidate value.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
EquationOfState<Dimension>::applyPressureLimits(Scalar P,
                                                const Scalar Di) const {
  REQUIRE(Di >= 0.0 and Di <= 1.0);

  // First apply standard limiting
  P = (P < mMinimumPressure ? (mMinPressureType == MaterialPressureMinType::PressureFloor ? 
                               mMinimumPressure :
                               0.0) :
          P > mMaximumPressure ? 
          mMaximumPressure :
          P);

  // Now any damage application to negative pressures
  if (P < 0.0) {
    P = (1.0 - Di)*P + Di*min(0.0, mMinimumDamagePressure);
  }

  return P;
}

}
