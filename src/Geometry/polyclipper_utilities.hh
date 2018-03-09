//------------------------------------------------------------------------------
// A collection of internal utilities common in PolyClipper.
//
// These largely repoduce things in Spheral to make PolyClipper standalone.
//------------------------------------------------------------------------------
#ifndef __polyclipper_utilities__
#define __polyclipper_utilities__

namespace PolyClipper {

//------------------------------------------------------------------------------
// Return the sign of the argument determined as follows:
//   
//    x > 0 -> sgn0(x) =  1
//    x = 0 -> sgn0(x) =  0
//    x < 0 -> sgn0(x) = -1
//------------------------------------------------------------------------------
inline
double
sgn0(const double x) {
  return (x > 0.0 ?  1.0 :
          x < 0.0 ? -1.0 :
          0.0);
}

//------------------------------------------------------------------------------
// Return the sign of the argument determined as follows:
//   
//    x >= 0 -> sgn(x) =  1
//    x <  0 -> sgn(x) = -1
//------------------------------------------------------------------------------
inline
double
sgn(const double x) {
  return (x >= 0.0 ?  1.0 : -1.0);
}

//------------------------------------------------------------------------------
// Fuzzy comparisons.
//------------------------------------------------------------------------------
template<typename DataType>
inline
bool fuzzyEqual(const DataType& lhs, const DataType& rhs,
                const double fuzz = 1.0e-15) {
  return std::abs(lhs - rhs)/std::max(1.0, std::abs(lhs) + std::abs(rhs)) < fuzz;
}

template<typename DataType>
inline
bool fuzzyLessThanOrEqual(const DataType& lhs, const DataType& rhs,
                          const double fuzz = 1.0e-15) {
  return lhs < rhs || fuzzyEqual(lhs, rhs, fuzz);
}

template<typename DataType>
inline
bool fuzzyGreaterThanOrEqual(const DataType& lhs, const DataType& rhs,
                             const double fuzz = 1.0e-15) {
  return lhs > rhs || fuzzyEqual(lhs, rhs, fuzz);
}

template<typename DataType>
inline
bool distinctlyLessThan(const DataType& lhs, const DataType& rhs,
                        const double fuzz = 1.0e-15) {
  return lhs < rhs && !fuzzyEqual(lhs, rhs, fuzz);
}

template<typename DataType>
inline
bool distinctlyGreaterThan(const DataType& lhs, const DataType& rhs,
                           const double fuzz = 1.0e-15) {
  return lhs > rhs && !fuzzyEqual(lhs, rhs, fuzz);
}

//------------------------------------------------------------------------------
// safeInv
//
// Take the inverse of a number, turning it to zero if the argument is 0.0.
//------------------------------------------------------------------------------
template<typename Value>
inline
Value
safeInv(const Value& x,
        const double fuzz = 1.0e-30) {
  return sgn(x)/std::max(fuzz, std::abs(x));
}

}

#endif
