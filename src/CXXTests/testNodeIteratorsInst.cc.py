text = """

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "CXXTests/testNodeIterators.cc"

namespace Spheral {
template string testGlobalAllNodeIterators<Dim< %(ndim)s > >(const DataBase<Dim< %(ndim)s > >&);

template string testGlobalInternalNodeIterators<Dim< %(ndim)s > >(const DataBase<Dim< %(ndim)s > >&);

template string testGlobalGhostNodeIterators<Dim< %(ndim)s > >(const DataBase<Dim< %(ndim)s > >&);

template string testGlobalMasterNodeIterators<Dim< %(ndim)s > >(const DataBase<Dim< %(ndim)s > >&);

template string testGlobalCoarseNodeIterators<Dim< %(ndim)s > >(const DataBase<Dim< %(ndim)s > >&);

template string testGlobalRefineNodeIterators<Dim< %(ndim)s > >(const DataBase<Dim< %(ndim)s > >&);
}

"""
