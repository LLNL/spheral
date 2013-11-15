//------------------------------------------------------------------------------
// Helper method for the GenerateCylindricalNodeDistribution3d node generator
// to generate the spun node distribution.
//------------------------------------------------------------------------------
#ifndef __Spheral_generateCylDistributionFromRZ__
#define __Spheral_generateCylDistributionFromRZ__

#include <vector>
#include "Geometry/Dimension.hh"

namespace Spheral {

void
generateCylDistributionFromRZ(std::vector<double>& x,
                              std::vector<double>& y,
                              std::vector<double>& z,
                              std::vector<double>& m,
                              std::vector<Dim<3>::SymTensor>& H,
                              std::vector<int>& globalIDs,
                              std::vector<std::vector<double> >& extraFields,
                              const double nNodePerh,
                              const double kernelExtent,
                              const double phi,
                              const int procID,
                              const int nProcs);

}


#endif
