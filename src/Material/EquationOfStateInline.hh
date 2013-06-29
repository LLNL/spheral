namespace Spheral {
namespace Material {

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

}
}
