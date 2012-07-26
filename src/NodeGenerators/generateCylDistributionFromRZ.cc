//------------------------------------------------------------------------------
// Helper method for the GenerateCylindricalNodeDistribution3d node generator
// to generate the spun node distribution.
//------------------------------------------------------------------------------
#ifdef USE_MPI
#include "mpi.h"
#endif

#include <vector>
#include <algorithm>

#include "Boundary/CylindricalBoundary.hh"
#include "Utilities/DBC.hh"
#include "Geometry/Dimension.hh"

using namespace std;
using namespace Spheral::BoundarySpace;

namespace Spheral {

void
generateCylDistributionFromRZ(vector<double>& x,
                              vector<double>& y,
                              vector<double>& z,
                              vector<double>& m,
                              vector<Dim<3>::SymTensor>& H,
                              vector<int>& globalIDs,
                              vector<vector<double> >& extraFields,
                              const double nNodePerh,
                              const double kernelExtent,
                              const double phi,
                              const int procID,
                              const int nProcs) {

  typedef Dim<3>::Vector Vector;
  typedef Dim<3>::Tensor Tensor;
  typedef Dim<3>::SymTensor SymTensor;

  // Pre-conditions.
  const int n = x.size();
  const int nextra = extraFields.size();
  VERIFY(y.size() == n and
         z.size() == n and
         m.size() == n and
         H.size() == n);
  for (int i = 0; i != nextra; ++i) VERIFY(extraFields[i].size() == n);
  VERIFY(count(z.begin(), z.end(), 0.0) == n);
  VERIFY(nNodePerh > 0.0);
  VERIFY(kernelExtent > 0.0);
  VERIFY(phi > 0.0);
  VERIFY(nProcs >= 1);
  VERIFY(procID >= 0 and procID < nProcs);

  // Make an initial pass to determine how many nodes we're going
  // to generate.
  int ntot = 0;
  for (int i = 0; i != n; ++i) {
    const SymTensor Hi = H[i];
    const double hzi = 1.0/(Hi*Vector(0.0, 0.0, 1.0)).magnitude();
    const double yi = y[i];
    double dphi = CylindricalBoundary::angularSpacing(yi, hzi, nNodePerh, kernelExtent);
    CHECK(distinctlyGreaterThan(dphi, 0.0));
    const int nsegment = max(1, int(phi/dphi + 0.5));
    ntot += nsegment;
  }

  // Determine how the global IDs should be partitioned between processors.
  const int ndomain0 = ntot/nProcs;
  const int remainder = ntot % nProcs;
  VERIFY(remainder < nProcs);
  const int ndomain = ndomain0 + (procID < remainder ? 1 : 0);
  const int minGlobalID = procID*ndomain0 + min(procID, remainder);
  const int maxGlobalID = minGlobalID + ndomain - 1;
  VERIFY(procID < nProcs - 1 || maxGlobalID == ntot - 1);
  
  // Copy the input.
  vector<double> xrz(x), yrz(y), zrz(z), mrz(m);
  vector<SymTensor> Hrz(H);
  vector<vector<double> > extrasrz(extraFields);

  // Prepare the lists we're going to rebuild.
  x = vector<double>();
  y = vector<double>();
  z = vector<double>();
  m = vector<double>();
  H = vector<SymTensor>();
  globalIDs = vector<int>();
  extraFields = vector<vector<double> >(nextra);

  // Iterate over the plane of input nodes, and rotate it out for the full 3-D 
  // distribution.
  int globalID = 0;
  for (int i = 0; i != n; ++i) {
    const SymTensor Hi = Hrz[i];
    const double hzi = 1.0/(Hi*Vector(0.0, 0.0, 1.0)).magnitude();
    const double xi = xrz[i];
    const double yi = yrz[i];
    const double mi = mrz[i];
    const int nhoopsegment = max(1, int(phi/CylindricalBoundary::angularSpacing(yi, hzi, nNodePerh, kernelExtent) + 0.5));
    const double dphi = phi/nhoopsegment;
    double phii = 0.0;
    const Vector posi = Vector(xi, yi, 0.0);
    for (int ihoop = 0; ihoop != nhoopsegment; ++ihoop) {
      const double phii = ihoop*dphi;
      if (globalID >= minGlobalID and globalID <= maxGlobalID) {
        const double xj = xi;
        const double yj = yi*cos(phii);
        const double zj = yi*sin(phii);
        x.push_back(xj);
        y.push_back(yj);
        z.push_back(zj);
        m.push_back(mi/nhoopsegment * phi/(2.0*M_PI));
        globalIDs.push_back(globalID);
        const Vector posj = Vector(xj, yj, zj);
        const Tensor R = CylindricalBoundary::reflectOperator(posi, posj);
        H.push_back((R*Hi*R).Symmetric());
        for (int ikey = 0; ikey != nextra; ++ikey) extraFields[ikey].push_back(extrasrz[ikey][i]);
      }
      ++globalID;
    }
  }

  // Post-conditions.
  VERIFY(x.size() == ndomain and
         y.size() == ndomain and
         z.size() == ndomain and
         m.size() == ndomain and
         globalIDs.size() == ndomain and
         H.size() == ndomain);
  for (int ikey = 0; ikey != nextra; ++ikey) VERIFY(extraFields[ikey].size() == ndomain);
#ifdef USE_MPI
  {
    int nlocal = x.size();
    int nglobal;
    MPI_Allreduce(&nlocal, &nglobal, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    VERIFY(nglobal == ntot);
  }
#endif

}

}
