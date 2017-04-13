#ifndef __PBGWRAPS_NODEGENERATORTYPES__
#define __PBGWRAPS_NODEGENERATORTYPES__

#include "NodeGenerators/generateCylDistributionFromRZ.hh"
#include "NodeGenerators/fillFacetedVolume.hh"
#include "NodeGenerators/relaxNodeDistribution.hh"
#include "NodeGenerators/readSiloPolyMesh.hh"
#include "NodeGenerators/centroidalRelaxNodesImpl.hh"
#include "NodeGenerators/compactFacetedVolumes.hh"
#include "NodeGenerators/chooseRandomNonoverlappingCenter.hh"
#include "Geometry/Dimension.hh"

namespace Spheral {

typedef WeightingFunctor<Dim<2> > WeightingFunctor2d;
typedef WeightingFunctor<Dim<3> > WeightingFunctor3d;

}

#endif
