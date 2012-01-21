#ifndef DataNameTraits_hh
#define DataNameTraits_hh

#include <string>
using namespace std;

namespace Spheral {
  template<typename DataType>
  struct DataNameTraits {
    static const std::string DataName;
  };
}

#else

namespace Spheral {
  template<typename DataType> struct DataNameTraits;
}

#endif
