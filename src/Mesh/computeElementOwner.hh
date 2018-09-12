//---------------------------------Spheral++----------------------------------//
// computeElementOwner
//
// Utility method to reduce the shared elements to the lowest common processor.
//
// Created by JMO, Mon Oct 29 22:01:30 PDT 2012
//----------------------------------------------------------------------------//
#ifndef __Spheral_computeElementOwner__
#define __Spheral_computeElementOwner__

#include <vector>

namespace Spheral {

std::vector<unsigned>
computeElementOwner(const unsigned numElements,
                    const std::vector<unsigned>& neighborDomains,
                    const std::vector<std::vector<unsigned> >& sharedElements);

}

#endif
