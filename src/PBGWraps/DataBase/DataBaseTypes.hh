#ifndef __PBGWRAPS_DATABASETYPES__
#define __PBGWRAPS_DATABASETYPES__

#include "Geometry/Dimension.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/StateBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Neighbor/ConnectivityMap.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Names!
//------------------------------------------------------------------------------
typedef StateBase<Dim<1> > StateBase1d;
typedef State<Dim<1> > State1d;
typedef StateDerivatives<Dim<1> > StateDerivatives1d;
typedef State<Dim<1> > State1d;

typedef StateBase<Dim<2> > StateBase2d;
typedef State<Dim<2> > State2d;
typedef StateDerivatives<Dim<2> > StateDerivatives2d;
typedef State<Dim<2> > State2d;

typedef StateBase<Dim<3> > StateBase3d;
typedef State<Dim<3> > State3d;
typedef StateDerivatives<Dim<3> > StateDerivatives3d;
typedef State<Dim<3> > State3d;

typedef DataBase<Dim<1> > DataBase1d;
typedef DataBase<Dim<2> > DataBase2d;
typedef DataBase<Dim<3> > DataBase3d;

//------------------------------------------------------------------------------
// Get fields as references from StateBase.
//------------------------------------------------------------------------------
template<typename Dimension, typename Value>
inline
Field<Dimension, Value>*
fieldFromStateBase(StateBase<Dimension>& self,
                   const typename StateBase<Dimension>::KeyType& key) {
  return &(self.field(key, Value()));
}

//------------------------------------------------------------------------------
// Get the connectivity map as a pointer from the DataBase.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
ConnectivityMap<Dimension>*
connectivityMapFromDataBase(const DataBase<Dimension>& db,
                            const bool buildGhostConnectivity) {
  return const_cast<ConnectivityMap<Dimension>*>(&(db.connectivityMap(buildGhostConnectivity)));
}

}

#endif
