#ifndef __PBGWRAPS_PHYSICSTYPES__
#define __PBGWRAPS_PHYSICSTYPES__

#include "Geometry/Dimension.hh"
#include "Physics/Physics.hh"
#include "Physics/GenericHydro.hh"
#include "Physics/GenericBodyForce.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Names!
//------------------------------------------------------------------------------
typedef Physics<Dim<1> > Physics1d;
typedef Physics<Dim<2> > Physics2d;
typedef Physics<Dim<3> > Physics3d;

typedef GenericHydro<Dim<1> > GenericHydro1d;
typedef GenericHydro<Dim<2> > GenericHydro2d;
typedef GenericHydro<Dim<3> > GenericHydro3d;

typedef GenericBodyForce<Dim<1> > GenericBodyForce1d;
typedef GenericBodyForce<Dim<2> > GenericBodyForce2d;
typedef GenericBodyForce<Dim<3> > GenericBodyForce3d;
    
}

typedef std::vector<Spheral::Physics<Spheral::Dim<1> >*> vector_of_Physics1d;
typedef std::vector<Spheral::Physics<Spheral::Dim<2> >*> vector_of_Physics2d;
typedef std::vector<Spheral::Physics<Spheral::Dim<3> >*> vector_of_Physics3d;

#endif
