//---------------------------------Spheral++----------------------------------//
// VolumeIntegrationFunctions.
// A set of helper methods to numerically evaluate the appropriate kernel
// normalization.
//
// Created by JMO, Mon Jan  6 16:24:19 PST 2003
//----------------------------------------------------------------------------//

namespace Spheral {

// Use Simpsons rule to evaluate the volume integral of the given kernel.
template<typename Dimension, typename KernelType>
double simpsonsVolumeIntegral(const KernelType& W,
                              const double rMin,
                              const double rMax,
                              const int numBins);

}

