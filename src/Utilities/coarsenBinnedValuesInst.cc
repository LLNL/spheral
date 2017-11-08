//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "coarsenBinnedValues.cc"

namespace Spheral {

  template void coarsenBinnedValues(vector<vector<double> >& values, const unsigned nxFine);
  // template void coarsenBinnedValues(vector<vector<Dim<1>::Vector> >& values, const unsigned nxFine);

  template void coarsenBinnedValues(vector<vector<double> >& values, const unsigned nxFine, const unsigned nyFine);
  // template void coarsenBinnedValues(vector<vector<Dim<2>::Vector> >& values, const unsigned nxFine, const unsigned nyFine);

  template void coarsenBinnedValues(vector<vector<double> >& values, const unsigned nxFine, const unsigned nyFine, const unsigned nzFine);
  // template void coarsenBinnedValues(vector<vector<Dim<3>::Vector> >& values, const unsigned nxFine, const unsigned nyFine, const unsigned nzFine);

}
