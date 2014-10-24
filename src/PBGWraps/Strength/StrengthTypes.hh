#ifndef __PBGWRAPS_STRENGTHTYPES__
#define __PBGWRAPS_STRENGTHTYPES__

#include "Geometry/Dimension.hh"
#include "Strength/SolidFieldNames.hh"
#include "Strength/SolidNodeList.hh"

//------------------------------------------------------------------------------
// Names!
//------------------------------------------------------------------------------
namespace Spheral {
namespace SolidMaterial {

typedef SolidNodeList<Dim<1> > SolidNodeList1d;
typedef SolidNodeList<Dim<2> > SolidNodeList2d;
typedef SolidNodeList<Dim<3> > SolidNodeList3d;

}
}

#endif
