//---------------------------------Spheral++----------------------------------//
// medianPosition
//
// Compute a definition of the median position for a collection of Vectors.
//
// Created by JMO, Thu Feb 18 11:25:55 PST 2010
//----------------------------------------------------------------------------//
#ifndef __Spheral_medianPosition__
#define __Spheral_medianPosition__

#include <vector>

namespace Spheral {

//------------------------------------------------------------------------------
// medianPosition
// Note this methd *will* shuffle the order of the input vector!
//------------------------------------------------------------------------------
template<typename Vector>
Vector
medianPosition(std::vector<Vector>& positions);

}

#endif
