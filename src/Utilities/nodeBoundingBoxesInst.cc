//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "nodeBoundingBoxes.cc"

namespace Spheral {

  template Field<Dim<1>, std::pair<Dim<1>::Vector, Dim<1>::Vector> > nodeBoundingBoxes(const NodeList<Dim<1> >& nodes);
  template Field<Dim<2>, std::pair<Dim<2>::Vector, Dim<2>::Vector> > nodeBoundingBoxes(const NodeList<Dim<2> >& nodes);
  template Field<Dim<3>, std::pair<Dim<3>::Vector, Dim<3>::Vector> > nodeBoundingBoxes(const NodeList<Dim<3> >& nodes);

  template FieldList<Dim<1>, std::pair<Dim<1>::Vector, Dim<1>::Vector> > nodeBoundingBoxes(const DataBase<Dim<1> >& dataBase);
  template FieldList<Dim<2>, std::pair<Dim<2>::Vector, Dim<2>::Vector> > nodeBoundingBoxes(const DataBase<Dim<2> >& dataBase);
  template FieldList<Dim<3>, std::pair<Dim<3>::Vector, Dim<3>::Vector> > nodeBoundingBoxes(const DataBase<Dim<3> >& dataBase);

}
