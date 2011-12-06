// BPL includes.
#include "boost/python.hpp"

#include "ExtendGenerateCylindricalDistribution3d.hh"

namespace Spheral {

BOOST_PYTHON_MODULE(NodeGenerators) {
  def("generateCylDistributionFromRZ", generateCylDistributionFromRZ);
}

}

