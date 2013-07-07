//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "orientedBoundingBox.cc"

namespace Spheral {

  template void orientedBoundingBox<Dim<1> >(const DataBase<Dim<1> >& dataBase, Dim<1>::Box& nodeBox, Dim<1>::Box& sampleBox);
  template void orientedBoundingBox<Dim<2> >(const DataBase<Dim<2> >& dataBase, Dim<2>::Box& nodeBox, Dim<2>::Box& sampleBox);
  template void orientedBoundingBox<Dim<3> >(const DataBase<Dim<3> >& dataBase, Dim<3>::Box& nodeBox, Dim<3>::Box& sampleBox);

}
