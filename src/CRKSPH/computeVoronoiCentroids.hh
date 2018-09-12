//---------------------------------Spheral++------------------------------------
// Compute centroids for each point using the Vornonoi tessellation.
//------------------------------------------------------------------------------
#ifndef __Spheral__computeVoronoiCentroids__
#define __Spheral__computeVoronoiCentroids__

#include "Geometry/Dimension.hh"
#include "Utilities/DBC.hh"
#include "Field/FieldList.hh"

namespace Spheral {

  inline
  FieldList<Dim<1>, Dim<1>::Vector>
  computeVoronoiCentroids(const FieldList<Dim<1>, Dim<1>::Vector>& position) { VERIFY2(false, "Unimplemented"); }

  FieldList<Dim<2>, Dim<2>::Vector>
  computeVoronoiCentroids(const FieldList<Dim<2>, Dim<2>::Vector>& position);

  inline
  FieldList<Dim<3>, Dim<3>::Vector>
  computeVoronoiCentroids(const FieldList<Dim<3>, Dim<3>::Vector>& position) { VERIFY2(false, "Unimplemented"); }

}

#endif
