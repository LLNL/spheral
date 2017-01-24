//------------------------------------------------------------------------------
// Helper method for the SiloPolyMeshGenerator which reads a polyhedral mesh 
// from a silo file and returns the geometry.
//------------------------------------------------------------------------------
#ifndef __Spheral_readSiloPolyMesh__
#define __Spheral_readSiloPolyMesh__

#include <vector>
#include "Geometry/Dimension.hh"

namespace Spheral {

void
readSiloPolyMesh(const std::string& fileName,
                 const std::string& meshName,
                 std::vector<Dim<3>::Vector>& positions,
                 std::vector<double>& volumes,
                 std::vector<Dim<3>::SymTensor>& H);

}


#endif
