//---------------------------------Spheral++----------------------------------//
// FieldBase -- An abstract base class to provide a generic handle on Fields.
//
// Created by JMO, Sat Oct 30 23:38:19 PDT 1999
//----------------------------------------------------------------------------//
#include "FieldBase.hh"
#include "FieldListBase.hh"
#include "Utilities/DBC.hh"

namespace Spheral {
namespace FieldSpace {

// //------------------------------------------------------------------------------
// // The clone method.
// // We actually wanted this to be a pure virtual method, but I'm not sure how
// // to handle it in python binding so punting for now.
// //------------------------------------------------------------------------------
// template<typename Dimension>
// boost::shared_ptr<FieldBase<Dimension> >
// FieldBase<Dimension>::
// clone() const {
//   VERIFY2(false, "FieldBase::clone ERROR -- consider this sucker pure virtual!");
// }

// //------------------------------------------------------------------------------
// // Trip the flag indicating there are new coarse nodes.
// //------------------------------------------------------------------------------
// template<typename Dimension>
// void
// FieldBase<Dimension>::
// notifyNewCoarseNodes() const {
//   for (std::vector<const FieldListBase*>::const_iterator itr = mFieldListBaseList.begin();
//        itr != mFieldListBaseList.end();
//        ++itr) (*itr)->notifyNewCoarseNodes();
// }

// //------------------------------------------------------------------------------
// // Trip the flag indicating there are new refine nodes.
// //------------------------------------------------------------------------------
// template<typename Dimension>
// void
// FieldBase<Dimension>::
// notifyNewRefineNodes() const {
//   for (std::vector<const FieldListBase*>::const_iterator itr = mFieldListBaseList.begin();
//        itr != mFieldListBaseList.end();
//        ++itr) (*itr)->notifyNewRefineNodes();
// }

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
template class FieldBase< Dim<1> >;
template class FieldBase< Dim<2> >;
template class FieldBase< Dim<3> >;
}
}
