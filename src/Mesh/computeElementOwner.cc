//---------------------------------Spheral++----------------------------------//
// computeElementOwner
//
// Utility method to reduce the shared elements to the lowest common processor.
//
// Created by JMO, Mon Oct 29 22:01:30 PDT 2012
//----------------------------------------------------------------------------//
#include "boost/foreach.hpp"
#include "Utilities/DBC.hh"

#include "computeElementOwner.hh"

namespace Spheral {
namespace MeshSpace {

using namespace std;

vector<unsigned>
computeElementOwner(const unsigned numElements,
                    const vector<unsigned>& neighborDomains,
                    const vector<vector<unsigned> >& sharedElements) {
  REQUIRE(neighborDomains.size() == sharedElements.size());

  const unsigned rank = Process::getRank();
  const unsigned numProcs = Process::getTotalNumberOfProcesses();

  // Prepare the result.
  vector<unsigned> result(numElements, rank);

  // There's really only work if we have more than one domain.
  if (numProcs > 1) {
    const unsigned numDomains = neighborDomains.size();
    for (unsigned idomain = 0; idomain != numDomains; ++idomain) {
      const unsigned neighbor = neighborDomains[idomain];
      BOOST_FOREACH(unsigned element, sharedElements[idomain]) {
        CHECK(element < numElements);
        result[element] = min(result[element], neighbor);
      }
    }
  }

  return result;
}

}
}
