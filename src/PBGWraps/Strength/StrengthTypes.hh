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

//------------------------------------------------------------------------------
// STL containers.
//------------------------------------------------------------------------------
typedef std::vector<Spheral::SolidMaterial::SolidNodeList<Dim<1> >*> vector_of_SolidNodeList1d;
typedef std::vector<Spheral::SolidMaterial::SolidNodeList<Dim<2> >*> vector_of_SolidNodeList2d;
typedef std::vector<Spheral::SolidMaterial::SolidNodeList<Dim<3> >*> vector_of_SolidNodeList3d;

typedef std::vector<const Spheral::SolidMaterial::SolidNodeList<Dim<1> >*> vector_of_const_SolidNodeList1d;
typedef std::vector<const Spheral::SolidMaterial::SolidNodeList<Dim<2> >*> vector_of_const_SolidNodeList2d;
typedef std::vector<const Spheral::SolidMaterial::SolidNodeList<Dim<3> >*> vector_of_const_SolidNodeList3d;

typedef std::vector<Spheral::SolidMaterial::SolidNodeList<Dim<1> >*>::iterator vector_of_SolidNodeList1d_iterator;
typedef std::vector<Spheral::SolidMaterial::SolidNodeList<Dim<2> >*>::iterator vector_of_SolidNodeList2d_iterator;
typedef std::vector<Spheral::SolidMaterial::SolidNodeList<Dim<3> >*>::iterator vector_of_SolidNodeList3d_iterator;

#endif
