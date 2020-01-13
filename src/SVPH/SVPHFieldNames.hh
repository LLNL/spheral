//---------------------------------Spheral++----------------------------------//
// SVPHFieldNames -- A collection of standard Field names for the SVPH 
// physics package.
//
// Created by JMO, Wed Jan  8 11:15:42 PST 2020
//----------------------------------------------------------------------------//
#ifndef _Spheral_SVPHFieldNames_
#define _Spheral_SVPHFieldNames_

#include <string>

namespace Spheral {

struct SVPHFieldNames {
  static const std::string A_SVPH;
  static const std::string B_SVPH;
  static const std::string gradB_SVPH;
};

}

#else

namespace Spheral {
  string SVPHFieldNames;
}

#endif
