//------------------------------------------------------------------------------
// Helper method for the GenerateCylindricalNodeDistribution3d node generator
// to generate the spun node distribution.
//------------------------------------------------------------------------------
#ifndef __Spheral_ExtendGenerateCylindricalDistribution3d__
#define __Spheral_ExtendGenerateCylindricalDistribution3d__

#ifdef USE_MPI
#include "mpi.h"
#endif

#ifndef __GCCXML__
#include "boost/python.hpp"
#include "boost/python/detail/api_placeholder.hpp"
#else
#include "../fakeboost.hh"
#endif

#include "Boundary/CylindricalBoundary.hh"
#include "Utilities/DBC.hh"
#include "Geometry/Dimension.hh"

namespace Spheral {

using namespace boost;
using namespace BoundarySpace;

inline
void
generateCylDistributionFromRZ(python::list x,
                              python::list y,
                              python::list z,
                              python::list m,
                              python::list H,
                              python::list globalIDs,
                              python::list extraFields,
                              python::object nNodePerhObj,
                              python::object kernelExtentObj,
                              python::object phiObj,
                              python::object procIDObj,
                              python::object nProcsObj) {

  typedef Dim<3>::Vector Vector;
  typedef Dim<3>::Tensor Tensor;
  typedef Dim<3>::SymTensor SymTensor;

  // Pre-conditions.
  const int n = python::len(x);
  const double nNodePerh = python::extract<double>(nNodePerhObj);
  const double kernelExtent = python::extract<double>(kernelExtentObj);
  const double phi = python::extract<double>(phiObj);
  const int procID = python::extract<int>(procIDObj);
  const int nProcs = python::extract<int>(nProcsObj);
  const int nextra = python::len(extraFields);
  VERIFY(python::len(y) == n &&
         python::len(z) == n &&
         python::len(m) == n &&
         python::len(H) == n);
  for (int i = 0; i != nextra; ++i) VERIFY(python::len(extraFields[i]) == n);
  VERIFY(z.count(0.0) == n);
  VERIFY(nNodePerh > 0.0);
  VERIFY(kernelExtent > 0.0);
  VERIFY(phi > 0.0);
  VERIFY(nProcs >= 1);
  VERIFY(procID >= 0 && procID < nProcs);

  // Make an initial pass to determine how many nodes we're going
  // to generate.
  int ntot = 0;
  for (int i = 0; i != n; ++i) {
    const SymTensor Hi = python::extract<SymTensor>(H[i]);
    const double hzi = 1.0/(Hi*Vector(0.0, 0.0, 1.0)).magnitude();
    const double yi = python::extract<double>(y[i]);
    const double dphi = CylindricalBoundary::angularSpacing(yi, hzi, nNodePerh, kernelExtent);
    CHECK(distinctlyGreaterThan(dphi, 0.0));
    double phii = 0.0;
    while (phii < phi - 0.5*dphi) {
      ++ntot;
      phii += dphi;
    }
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
  python::list xrz, yrz, zrz, mrz, Hrz;
  xrz.extend(x);
  yrz.extend(y);
  zrz.extend(z);
  mrz.extend(m);
  Hrz.extend(H);
  python::list extrasrz;
  for (int i = 0; i != nextra; ++i) {
    extrasrz.append(python::list());
    VERIFY(python::len(extrasrz) == i + 1);
    python::list values = python::extract<python::list>(extraFields[i]);
    python::list valuesrz = python::extract<python::list>(extraFields[i]);
    valuesrz.extend(values);
    VERIFY(python::len(valuesrz) == n);
  }

  // Prepare the lists we're going to rebuild.
  for (int i = 0; i != n; ++i) {
    x.pop();
    y.pop();
    z.pop();
    m.pop();
    H.pop();
    globalIDs.pop();
  }
  VERIFY(python::len(x) == 0);
  VERIFY(python::len(y) == 0);
  VERIFY(python::len(z) == 0);
  VERIFY(python::len(m) == 0);
  VERIFY(python::len(H) == 0);
  VERIFY(python::len(globalIDs) == 0);
  for (int ikey = 0; ikey != nextra; ++ikey) {
    python::list values = python::extract<python::list>(extraFields[ikey]);
    for (int i = 0; i != n; ++i) values.pop();
    VERIFY(python::len(python::extract<python::list>(extraFields[ikey])) == 0);
  }

  // Iterate over the plane of input nodes, and rotate it out for the full 3-D 
  // distribution.
  int globalID = 0;
  for (int i = 0; i != n; ++i) {
    const SymTensor Hi = python::extract<SymTensor>(Hrz[i]);
    const double hzi = 1.0/(Hi*Vector(0.0, 0.0, 1.0)).magnitude();
    const double xi = python::extract<double>(xrz[i]);
    const double yi = python::extract<double>(yrz[i]);
    const double mi = python::extract<double>(mrz[i]);
    const int nhoopsegment = int(phi/CylindricalBoundary::angularSpacing(yi, hzi, nNodePerh, kernelExtent) + 0.5);
    const double dphi = phi/nhoopsegment;
    double phii = 0.0;
    const Vector posi = Vector(xi, yi, 0.0);
    for (int ihoop = 0; ihoop != nhoopsegment; ++ihoop) {
      const double phii = ihoop*dphi;
      if (globalID >= minGlobalID && globalID <= maxGlobalID) {
        const double xj = xi;
        const double yj = yi*cos(phii);
        const double zj = yi*sin(phii);
        x.append(xj);
        y.append(yj);
        z.append(zj);
        m.append(mi/nhoopsegment * phi/(2.0*M_PI));
        globalIDs.append(globalID);
        const Vector posj = Vector(xj, yj, zj);
        const Tensor R = CylindricalBoundary::reflectOperator(posi, posj);
        H.append((R*Hi*R).Symmetric());

        for (int ikey = 0; ikey != nextra; ++ikey) {
          python::list values = python::extract<python::list>(extraFields[ikey]);
          python::list valuesrz = python::extract<python::list>(extrasrz[ikey]);
          values.append(valuesrz[i]);
        }
      }
      ++globalID;
    }
  }

  // Post-conditions.
  VERIFY(python::len(x) == ndomain);
  VERIFY(python::len(y) == python::len(x) &&
         python::len(z) == python::len(x) &&
         python::len(m) == python::len(x) &&
         python::len(H) == python::len(x));
  for (int ikey = 0; ikey != nextra; ++ikey) 
    VERIFY(python::len(extraFields[ikey]) == python::len(x));
#ifdef USE_MPI
  {
    int nlocal = python::len(x);
    int nglobal;
    MPI_Allreduce(&nlocal, &nglobal, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    VERIFY(nglobal == ntot);
  }
#endif

}

}

#endif
