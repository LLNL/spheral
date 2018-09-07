#ifndef __PBGWRAPS_PARMETISTYPES__
#define __PBGWRAPS_PARMETISTYPES__

#include "Distributed/ParmetisRedistributeNodes.hh"

//------------------------------------------------------------------------------
// Names!
//------------------------------------------------------------------------------
namespace Spheral {

typedef ParmetisRedistributeNodes<Dim<1> > ParmetisRedistributeNodes1d;
typedef ParmetisRedistributeNodes<Dim<2> > ParmetisRedistributeNodes2d;
typedef ParmetisRedistributeNodes<Dim<3> > ParmetisRedistributeNodes3d;

}

#endif
