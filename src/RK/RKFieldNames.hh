//---------------------------------Spheral++----------------------------------//
// RKFieldNames -- A collection of standard Field names for the Reproducing
// Kernel (RK) package.
//
// Created by JMO, Sat Dec 28 13:18:36 PST 2019
//----------------------------------------------------------------------------//
#ifndef _Spheral_RKFieldNames_
#define _Spheral_RKFieldNames_

#include "RK/RKCorrectionParams.hh"
#include <string>

namespace Spheral {

struct RKFieldNames {
  static const std::string rkOrders;
  static const std::string rkCorrections(const RKOrder order);
  static const std::string reproducingKernel(const RKOrder order);
};

}

#endif
