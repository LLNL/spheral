//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "spheralWildMagicConverters.cc"

namespace Spheral {

  template vector<Dim<1>::WMVector> wildMagicPositions<Dim<1> >(const NodeList<Dim<1> >& nodes);
  template vector<Dim<2>::WMVector> wildMagicPositions<Dim<2> >(const NodeList<Dim<2> >& nodes);
  template vector<Dim<3>::WMVector> wildMagicPositions<Dim<3> >(const NodeList<Dim<3> >& nodes);

  template vector<Dim<1>::WMVector> wildMagicPositions<Dim<1> >(const DataBase<Dim<1> >& dataBase);
  template vector<Dim<2>::WMVector> wildMagicPositions<Dim<2> >(const DataBase<Dim<2> >& dataBase);
  template vector<Dim<3>::WMVector> wildMagicPositions<Dim<3> >(const DataBase<Dim<3> >& dataBase);

}
