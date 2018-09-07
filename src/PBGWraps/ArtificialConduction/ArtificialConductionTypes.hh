#ifndef __PBGWRAPS_ARTIFICIALCONDUCTIONTYPES__
#define __PBGWRAPS_ARTIFICIALCONDUCTIONTYPES__

#include "Geometry/Dimension.hh"
#include "ArtificialConduction/ArtificialConduction.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Names!
//------------------------------------------------------------------------------

typedef ArtificialConduction<Dim<1> > ArtificialConduction1d;
typedef ArtificialConduction<Dim<2> > ArtificialConduction2d;
typedef ArtificialConduction<Dim<3> > ArtificialConduction3d;
    
}

#endif
