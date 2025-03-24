//---------------------------------Spheral++----------------------------------//
// SpheralMessage
//
// Utilities to facilitate messages and warnings.
//
// Created by JMO, Mon Dec 16 15:20:42 PST 2024
//----------------------------------------------------------------------------//
#ifndef __Spheral_SpheralMessage__
#define __Spheral_SpheralMessage__

#include "Distributed/Process.hh"
#include <string>
#include <iostream>

namespace Spheral {

inline void SpheralMessage(const std::string& msg) {
  if (Process::getRank() == 0) std::cerr << msg << std::endl;
}

inline void DeprecationWarning(const std::string& msg) {
  SpheralMessage("DEPRECATION Warning: " + msg);
}

}

#endif

