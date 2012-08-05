//---------------------------------Spheral++----------------------------------//
// SpheralError -- a collection of routines to help handle error conditions
//                 in Spheral.
//
// Created by JMO, Tue Jan 11 22:49:52 PST 2000
//----------------------------------------------------------------------------//

#include "SpheralError.hh"

#include <iostream>
using namespace std;

namespace Spheral {

void SpheralError(const string& errorMessage) {
  throw SpheralWarning(errorMessage);
}

}
