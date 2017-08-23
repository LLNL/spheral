text = """

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "testNodeIterators.cc"

namespace Spheral {
  namespace Testing {
    template string testGlobalAllNodeIterators<Dim< %(ndim)s > >(const DataBase<Dim< %(ndim)s > >&);

    template string testGlobalInternalNodeIterators<Dim< %(ndim)s > >(const DataBase<Dim< %(ndim)s > >&);

    template string testGlobalGhostNodeIterators<Dim< %(ndim)s > >(const DataBase<Dim< %(ndim)s > >&);

    template string testGlobalMasterNodeIterators<Dim< %(ndim)s > >(const DataBase<Dim< %(ndim)s > >&);

    template string testGlobalCoarseNodeIterators<Dim< %(ndim)s > >(const DataBase<Dim< %(ndim)s > >&);

    template string testGlobalRefineNodeIterators<Dim< %(ndim)s > >(const DataBase<Dim< %(ndim)s > >&);
  }
}

"""
