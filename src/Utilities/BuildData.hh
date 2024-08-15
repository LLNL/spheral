#ifndef SPHERAL_BUILD_DATA_HPP
#define SPHERAL_BUILD_DATA_HPP

#include "Spheral/config.hh"
#include <string>

namespace Spheral {

class BuildData {
  public:
    static const std::string cxx_compiler_id;
    //static const std::string get_cxx_compiler_id(){ return cxx_compiler_id; };
};

} //  namespace Spheral

#endif //  SPHERAL_BUILD_DATA_HPP

