#ifndef PRINTGENERICOBJECT_HH
#define PRINTGENERICOBJECT_HH
//------------------------------------------------------------------------------
// Generic pattern for a string representation of a Spheral++ object.
//------------------------------------------------------------------------------
#include <string>
#include <sstream>

namespace Spheral {

template<typename ObjectType>
std::string
printObject(const ObjectType& self, const std::string& label) {
  std::stringstream result;
  result << label << "<" << &self << ">";
  return result.str();
}

}

#else
namespace Spheral {
  template<tyepname ObjectType>
  std::string printObject(const ObjectType& self, const std::string& label);
}
#endif
