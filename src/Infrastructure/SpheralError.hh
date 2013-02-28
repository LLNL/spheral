//---------------------------------Spheral++----------------------------------//
// SpheralError -- a collection of routines to help handle error conditions
//                 in Spheral.
//
// Created by JMO, Tue Jan 11 22:49:52 PST 2000
//----------------------------------------------------------------------------//

#ifndef __Spheral_SpheralError_hh__
#define __Spheral_SpheralError_hh__

#include <iostream>
#include <string>

namespace Spheral {

class SpheralError {
public:
  SpheralError():
    mErrorMessage("General Spheral error condition.") { printError(); }
  SpheralError(std::string errorMessage):
    mErrorMessage(errorMessage) { printError(); }

  void printError() {
    std::cerr << mErrorMessage << std::endl;
  }

private:
  std::string mErrorMessage;
};

}

#else

// Forward declaration.
namespace Spheral {
  class SpheralError;
}

#endif
