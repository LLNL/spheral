#ifndef __PBGWRAPS_ARTIFICIALVISCOSITYTYPES__
#define __PBGWRAPS_ARTIFICIALVISCOSITYTYPES__

#include "Geometry/Dimension.hh"
#include "ArtificialViscosity/ArtificialViscosity.hh"
#include "ArtificialViscosity/MonaghanGingoldViscosity.hh"
#include "ArtificialViscosity/CRKSPHMonaghanGingoldViscosity.hh"
#include "ArtificialViscosity/MorrisMonaghanReducingViscosity.hh"
#include "ArtificialViscosity/CullenDehnenViscosity.hh"
#include "ArtificialViscosity/TensorMonaghanGingoldViscosity.hh"
#include "ArtificialViscosity/FiniteVolumeViscosity.hh"
#include "ArtificialViscosity/TensorSVPHViscosity.hh"
#include "ArtificialViscosity/TensorCRKSPHViscosity.hh"
#include "ArtificialViscosity/VonNeumanViscosity.hh"
#include "ArtificialViscosity/MonaghanGingoldViscosityGSRZ.hh"
//#include "ArtificialViscosity/MonaghanGingoldViscosityRZ.hh"
//#include "ArtificialViscosity/CRKSPHMonaghanGingoldViscosityRZ.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Names!
//------------------------------------------------------------------------------
typedef ArtificialViscosity<Dim<1> > ArtificialViscosity1d;
typedef ArtificialViscosity<Dim<2> > ArtificialViscosity2d;
typedef ArtificialViscosity<Dim<3> > ArtificialViscosity3d;

typedef MonaghanGingoldViscosity<Dim<1> > MonaghanGingoldViscosity1d;
typedef MonaghanGingoldViscosity<Dim<2> > MonaghanGingoldViscosity2d;
typedef MonaghanGingoldViscosity<Dim<3> > MonaghanGingoldViscosity3d;
    
typedef CRKSPHMonaghanGingoldViscosity<Dim<1> > CRKSPHMonaghanGingoldViscosity1d;
typedef CRKSPHMonaghanGingoldViscosity<Dim<2> > CRKSPHMonaghanGingoldViscosity2d;
typedef CRKSPHMonaghanGingoldViscosity<Dim<3> > CRKSPHMonaghanGingoldViscosity3d;
    
typedef MorrisMonaghanReducingViscosity<Dim<1> > MorrisMonaghanReducingViscosity1d;
typedef MorrisMonaghanReducingViscosity<Dim<2> > MorrisMonaghanReducingViscosity2d;
typedef MorrisMonaghanReducingViscosity<Dim<3> > MorrisMonaghanReducingViscosity3d;

typedef CullenDehnenViscosity<Dim<1> > CullenDehnenViscosity1d;
typedef CullenDehnenViscosity<Dim<2> > CullenDehnenViscosity2d;
typedef CullenDehnenViscosity<Dim<3> > CullenDehnenViscosity3d;

typedef TensorMonaghanGingoldViscosity<Dim<1> > TensorMonaghanGingoldViscosity1d;
typedef TensorMonaghanGingoldViscosity<Dim<2> > TensorMonaghanGingoldViscosity2d;
typedef TensorMonaghanGingoldViscosity<Dim<3> > TensorMonaghanGingoldViscosity3d;

typedef FiniteVolumeViscosity<Dim<1> > FiniteVolumeViscosity1d;
typedef FiniteVolumeViscosity<Dim<2> > FiniteVolumeViscosity2d;
typedef FiniteVolumeViscosity<Dim<3> > FiniteVolumeViscosity3d;

typedef TensorSVPHViscosity<Dim<1> > TensorSVPHViscosity1d;
typedef TensorSVPHViscosity<Dim<2> > TensorSVPHViscosity2d;
typedef TensorSVPHViscosity<Dim<3> > TensorSVPHViscosity3d;

typedef TensorCRKSPHViscosity<Dim<1> > TensorCRKSPHViscosity1d;
typedef TensorCRKSPHViscosity<Dim<2> > TensorCRKSPHViscosity2d;
typedef TensorCRKSPHViscosity<Dim<3> > TensorCRKSPHViscosity3d;

typedef VonNeumanViscosity<Dim<1> > VonNeumanViscosity1d;
typedef VonNeumanViscosity<Dim<2> > VonNeumanViscosity2d;
typedef VonNeumanViscosity<Dim<3> > VonNeumanViscosity3d;

}

#endif
