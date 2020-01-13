#include "RK/RKFieldNames.hh"

namespace Spheral {

const std::string RKFieldNames::rkOrders = "rkOrders";
const std::string RKFieldNames::rkCorrectionsBase = "rkCorrections_";
const std::string RKFieldNames::reproducingKernelBase = "reproducingKernel_";
const std::string RKFieldNames::rkCorrections(const RKOrder order)     { return RKFieldNames::rkCorrectionsBase + std::to_string(static_cast<int>(order)); }
const std::string RKFieldNames::reproducingKernel(const RKOrder order) { return RKFieldNames::reproducingKernelBase + std::to_string(static_cast<int>(order)); }
const RKOrder RKFieldNames::correctionOrder(const std::string& x)      { const auto i = x.find("_"); return static_cast<RKOrder>(std::stoi(x.substr(i+1))); }

}
