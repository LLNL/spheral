#include <iostream>
#include <iomanip>
#include <sstream>

#include "vectorstringUtilities.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// vector<string> -> string specialization.
//------------------------------------------------------------------------------
template<>
std::string
vector2string<std::string>(const std::vector<std::string>& val,
                           const int precision) {
  std::ostringstream ss;
  ss.precision(precision);
  ss << val.size() << '\0';
  for (std::vector<std::string>::const_iterator itr = val.begin();
       itr != val.end();
       ++itr) ss << *itr << '\0';
  ss << std::ends;
  return ss.str();
}

//------------------------------------------------------------------------------
// string -> vector<string> specialization.
//------------------------------------------------------------------------------
template<>
std::vector<std::string>
string2vector<std::string>(const std::string& val) {
  std::istringstream ss(val);
  std::string x;
  std::getline(ss, x, '\0');
  std::istringstream xs(x);
  int size;
  xs >> size;
  std::vector<std::string> result;
  result.reserve(size);
  while (std::getline(ss, x, '\0')) result.push_back(x);
  result.pop_back();
//   VERIFY(result.size() == size);
  return result;
}


}
