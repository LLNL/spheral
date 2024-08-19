//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Utilities/nodeBoundingBoxes.cc"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)

  template Field<Dim<1>, std::pair<Dim<1>::Vector, Dim<1>::Vector> > nodeBoundingBoxes(const NodeList<Dim<1> >& nodes);
  template FieldList<Dim<1>, std::pair<Dim<1>::Vector, Dim<1>::Vector> > nodeBoundingBoxes(const DataBase<Dim<1> >& dataBase);

#endif

#if defined(SPHERAL_ENABLE_2D)

  template Field<Dim<2>, std::pair<Dim<2>::Vector, Dim<2>::Vector> > nodeBoundingBoxes(const NodeList<Dim<2> >& nodes);
  template FieldList<Dim<2>, std::pair<Dim<2>::Vector, Dim<2>::Vector> > nodeBoundingBoxes(const DataBase<Dim<2> >& dataBase);

#endif

#if defined(SPHERAL_ENABLE_3D)

  template Field<Dim<3>, std::pair<Dim<3>::Vector, Dim<3>::Vector> > nodeBoundingBoxes(const NodeList<Dim<3> >& nodes);
  template FieldList<Dim<3>, std::pair<Dim<3>::Vector, Dim<3>::Vector> > nodeBoundingBoxes(const DataBase<Dim<3> >& dataBase);

#endif
}