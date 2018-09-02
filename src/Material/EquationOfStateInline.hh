namespace Spheral {

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
EquationOfState<Dimension>::applyPressureLimits(const Scalar P) const {
  return (P < mMinimumPressure ? (mMinPressureType == MaterialPressureMinType::PressureFloor ? 
                                  mMinimumPressure :
                                  0.0) :
          P > mMaximumPressure ? 
          mMaximumPressure :
          P);
}

}
