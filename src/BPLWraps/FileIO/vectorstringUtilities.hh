#include <string>
#include <sstream>
#include <vector>
#include "Geometry/Dimension.hh"

#include "DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Generic vector -> string conversion.
//------------------------------------------------------------------------------
template<typename T>
inline
std::string
vector2string(const std::vector<T>& val,
              const int precision) {
  std::ostringstream ss;
  ss.precision(precision);
  ss << val.size() << " ";
  for (typename std::vector<T>::const_iterator itr = val.begin();
       itr != val.end();
       ++itr) ss << *itr << " ";
  ss << std::ends;
  return ss.str();
}

//------------------------------------------------------------------------------
// Generic string -> vector conversion.
//------------------------------------------------------------------------------
template<typename T>
inline
std::vector<T>
string2vector(const std::string& val) {
  std::istringstream ss(val);
  int size;
  ss >> size;
  std::vector<T> result;
  result.reserve(size);
  T x;
  while (ss >> x) result.push_back(x);
  VERIFY(result.size() == size);
  return result;
}

//------------------------------------------------------------------------------
// Declare the specializations.
//------------------------------------------------------------------------------
template<> std::string vector2string(const std::vector<int>& val, const int precision);
template<> std::string vector2string(const std::vector<bool>& val, const int precision);
template<> std::string vector2string(const std::vector<double>& val, const int precision);
template<> std::string vector2string(const std::vector<std::string>& val, const int precision);
template<> std::string vector2string(const std::vector<Spheral::Dim<1>::Vector>& val, const int precision);
template<> std::string vector2string(const std::vector<Spheral::Dim<2>::Vector>& val, const int precision);
template<> std::string vector2string(const std::vector<Spheral::Dim<3>::Vector>& val, const int precision);
template<> std::string vector2string(const std::vector<Spheral::Dim<1>::Tensor>& val, const int precision);
template<> std::string vector2string(const std::vector<Spheral::Dim<2>::Tensor>& val, const int precision);
template<> std::string vector2string(const std::vector<Spheral::Dim<3>::Tensor>& val, const int precision);
template<> std::string vector2string(const std::vector<Spheral::Dim<1>::SymTensor>& val, const int precision);
template<> std::string vector2string(const std::vector<Spheral::Dim<2>::SymTensor>& val, const int precision);
template<> std::string vector2string(const std::vector<Spheral::Dim<3>::SymTensor>& val, const int precision);
template<> std::string vector2string(const std::vector<Spheral::Dim<1>::ThirdRankTensor>& val, const int precision);
template<> std::string vector2string(const std::vector<Spheral::Dim<2>::ThirdRankTensor>& val, const int precision);
template<> std::string vector2string(const std::vector<Spheral::Dim<3>::ThirdRankTensor>& val, const int precision);

template<> std::vector<int> string2vector<int>(const std::string& val);
template<> std::vector<bool> string2vector<bool>(const std::string& val);
template<> std::vector<double> string2vector<double>(const std::string& val);
template<> std::vector<std::string> string2vector<std::string>(const std::string& val);
template<> std::vector<Spheral::Dim<1>::Vector> string2vector<Spheral::Dim<1>::Vector>(const std::string& val);
template<> std::vector<Spheral::Dim<2>::Vector> string2vector<Spheral::Dim<2>::Vector>(const std::string& val);
template<> std::vector<Spheral::Dim<3>::Vector> string2vector<Spheral::Dim<3>::Vector>(const std::string& val);
template<> std::vector<Spheral::Dim<1>::Tensor> string2vector<Spheral::Dim<1>::Tensor>(const std::string& val);
template<> std::vector<Spheral::Dim<2>::Tensor> string2vector<Spheral::Dim<2>::Tensor>(const std::string& val);
template<> std::vector<Spheral::Dim<3>::Tensor> string2vector<Spheral::Dim<3>::Tensor>(const std::string& val);
template<> std::vector<Spheral::Dim<1>::SymTensor> string2vector<Spheral::Dim<1>::SymTensor>(const std::string& val);
template<> std::vector<Spheral::Dim<2>::SymTensor> string2vector<Spheral::Dim<2>::SymTensor>(const std::string& val);
template<> std::vector<Spheral::Dim<3>::SymTensor> string2vector<Spheral::Dim<3>::SymTensor>(const std::string& val);
template<> std::vector<Spheral::Dim<1>::ThirdRankTensor> string2vector<Spheral::Dim<1>::ThirdRankTensor>(const std::string& val);
template<> std::vector<Spheral::Dim<2>::ThirdRankTensor> string2vector<Spheral::Dim<2>::ThirdRankTensor>(const std::string& val);
template<> std::vector<Spheral::Dim<3>::ThirdRankTensor> string2vector<Spheral::Dim<3>::ThirdRankTensor>(const std::string& val);

}
