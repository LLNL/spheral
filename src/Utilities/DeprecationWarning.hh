//---------------------------------Spheral++----------------------------------//
// DeprecationWarning.
//
// Print a warning about a deprecated feature.
//
// Created by JMO, Mon Dec 16 15:20:42 PST 2024
//----------------------------------------------------------------------------//
#ifndef __Spheral_DeprecationWarning__
#define __Spheral_DeprecationWarning__

#include "Distributed/Process.hh"
#include <string>
#include <iostream>

namespace Spheral {
inline
void DeprecationWarning(const std::string& msg) {
  if (Process::getRank() == 0) std::cerr << "DEPRECATION Warning: " << msg << std::endl;
}
}

#endif

