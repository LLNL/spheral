//------------------------------------------------------------------------------
// Construct a flaw distribution for a set of nodes according to a Weibull 
// (power-law) distribution.
//------------------------------------------------------------------------------
#ifndef __Spheral_weibullFlawDistribution__
#define __Spheral_weibullFlawDistribution__

#ifndef __GCCXML__
#include <vector>
#else
#include "fakestl.hh"
#endif

// Foward declarations.
namespace Spheral {
  namespace NodeSpace {
    template<typename Dimension> class FluidNodeList;
  }
  namespace FieldSpace {
    template<typename Dimension, typename DataType> class Field;
  }
}

namespace Spheral {
namespace PhysicsSpace {

//------------------------------------------------------------------------------
// Implements the Benz-Asphaug algorithm, starting with a minimum based
// on the volume of the simulation.
//------------------------------------------------------------------------------
template<typename Dimension>
FieldSpace::Field<Dimension, std::vector<double> >
weibullFlawDistributionBenzAsphaug(double volume,
                                   const double volumeStretchFactor,
                                   const unsigned seed,
                                   const double kWeibull,
                                   const double mWeibull,
                                   const NodeSpace::FluidNodeList<Dimension>& nodeList,
                                   const int minFlawsPerNode,
                                   const int minTotalFlaws);

//------------------------------------------------------------------------------
// Implements the Owen algorithm, stochastically seeding flaws with a maximum
// value per node chosen based on the volume of the node.
//------------------------------------------------------------------------------
template<typename Dimension>
FieldSpace::Field<Dimension, std::vector<double> >
weibullFlawDistributionOwen(const unsigned seed,
                            const double kWeibull,
                            const double mWeibull,
                            const NodeSpace::FluidNodeList<Dimension>& nodeList,
                            const int numFlawsPerNode,
                            const double volumeMultiplier = 1.0);

}
}

#endif
