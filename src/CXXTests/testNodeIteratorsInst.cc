//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "testNodeIterators.cc"

namespace Spheral {
  namespace Testing {
    template string testGlobalAllNodeIterators<Dim<1> >(const DataBase<Dim<1> >&);
    template string testGlobalAllNodeIterators<Dim<2> >(const DataBase<Dim<2> >&);
    template string testGlobalAllNodeIterators<Dim<3> >(const DataBase<Dim<3> >&);

    template string testGlobalInternalNodeIterators<Dim<1> >(const DataBase<Dim<1> >&);
    template string testGlobalInternalNodeIterators<Dim<2> >(const DataBase<Dim<2> >&);
    template string testGlobalInternalNodeIterators<Dim<3> >(const DataBase<Dim<3> >&);

    template string testGlobalGhostNodeIterators<Dim<1> >(const DataBase<Dim<1> >&);
    template string testGlobalGhostNodeIterators<Dim<2> >(const DataBase<Dim<2> >&);
    template string testGlobalGhostNodeIterators<Dim<3> >(const DataBase<Dim<3> >&);

    template string testGlobalMasterNodeIterators<Dim<1> >(const DataBase<Dim<1> >&);
    template string testGlobalMasterNodeIterators<Dim<2> >(const DataBase<Dim<2> >&);
    template string testGlobalMasterNodeIterators<Dim<3> >(const DataBase<Dim<3> >&);

    template string testGlobalCoarseNodeIterators<Dim<1> >(const DataBase<Dim<1> >&);
    template string testGlobalCoarseNodeIterators<Dim<2> >(const DataBase<Dim<2> >&);
    template string testGlobalCoarseNodeIterators<Dim<3> >(const DataBase<Dim<3> >&);

    template string testGlobalRefineNodeIterators<Dim<1> >(const DataBase<Dim<1> >&);
    template string testGlobalRefineNodeIterators<Dim<2> >(const DataBase<Dim<2> >&);
    template string testGlobalRefineNodeIterators<Dim<3> >(const DataBase<Dim<3> >&);

  }
}
