#include "RK/RKFieldNames.hh"

namespace Spheral {

const std::string rkOrders = "rkOrders";
static const std::string rkCorrections(const RKOrder order)     { return "rkCorrections_" + std::to_string(static_cast<int>(order)); }
static const std::string reproducingKernel(const RKOrder order) { return "reproducingKernel_" + std::to_string(static_cast<int>(order)); }

}
