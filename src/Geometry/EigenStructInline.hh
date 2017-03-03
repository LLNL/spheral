#include <string>

namespace Spheral {

//******************************************************************************
// Global functions.
//******************************************************************************
//------------------------------------------------------------------------------
// Input (istream) operator.
//------------------------------------------------------------------------------
template<int nDim>
inline
std::istream&
operator>>(std::istream& is, Spheral::EigenStruct<nDim>& eigen) {
  std::string parenthesis;
  is >> parenthesis;
  is >> eigen.eigenValues;
  is >> eigen.eigenVectors;
  is >> parenthesis;
  return is;
}

//------------------------------------------------------------------------------
// Output (ostream) operator.
//------------------------------------------------------------------------------
template<int nDim>
inline
std::ostream&
operator<<(std::ostream& os, const Spheral::EigenStruct<nDim>& eigen) {
  os << "( " << eigen.eigenValues << " " << eigen.eigenVectors << " )";
  return os;
}

}
