//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "nodeBoundingBoxes.cc"

namespace Spheral {

  template Field<Dim<1>, Dim<1>::Box> nodeBoundingBoxes(const NodeList<Dim<1> >& nodes);
  template Field<Dim<2>, Dim<2>::Box> nodeBoundingBoxes(const NodeList<Dim<2> >& nodes);
  template Field<Dim<3>, Dim<3>::Box> nodeBoundingBoxes(const NodeList<Dim<3> >& nodes);

  template FieldList<Dim<1>, Dim<1>::Box> nodeBoundingBoxes(const DataBase<Dim<1> >& dataBase);
  template FieldList<Dim<2>, Dim<2>::Box> nodeBoundingBoxes(const DataBase<Dim<2> >& dataBase);
  template FieldList<Dim<3>, Dim<3>::Box> nodeBoundingBoxes(const DataBase<Dim<3> >& dataBase);

}
