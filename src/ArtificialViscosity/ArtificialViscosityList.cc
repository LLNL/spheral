//---------------------------------Spheral++----------------------------------//
// ArtificialViscosity -- Make a list of Q's which will sum to act like a 
// single Q.  When one viscosity just won't do!
//----------------------------------------------------------------------------//
#include "ArtificialViscosityList.hh"
#include "Kernel/TableKernel.hh"

#include "DBC.hh"
#include "cdebug.hh"

#include <algorithm>

namespace Spheral {
namespace ArtificialViscositySpace {
using namespace std;

//------------------------------------------------------------------------------
// Default constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
ArtificialViscosityList<Dimension>::
ArtificialViscosityList():
  ArtificialViscosity<Dimension>(),
  vector<ArtificialViscosity<Dimension>*>() {
  cdebug << "ArtificialViscosityList::ArtificialViscosityList()" << endl;
}

//------------------------------------------------------------------------------
// Construct with the given set of artificial viscosities
//------------------------------------------------------------------------------
template<typename Dimension>
ArtificialViscosityList<Dimension>::
ArtificialViscosityList(const vector<ArtificialViscosity<Dimension>*>& QPtrs):
  ArtificialViscosity<Dimension>(),
  vector<ArtificialViscosity<Dimension>*>(QPtrs) {
  cdebug << "ArtificialViscosityList::ArtificialViscosityList(vector)" << endl;
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
ArtificialViscosityList<Dimension>::
~ArtificialViscosityList() {
  cdebug << "ArtificialViscosityList::~ArtificialViscosityList()" << endl;
}

//------------------------------------------------------------------------------
// Initialize for the FluidNodeLists in the given DataBase.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ArtificialViscosityList<Dimension>::
initialize(const DataBase<Dimension>& dataBase,
           typename ArtificialViscosity<Dimension>::ConstBoundaryIterator boundaryBegin,
           typename ArtificialViscosity<Dimension>::ConstBoundaryIterator boundaryEnd,
	   const typename Dimension::Scalar time,
	   const typename Dimension::Scalar dt,
           const TableKernel<Dimension>& W) {
  cdebug << "ArtificialViscosityList::initialize()" << endl;

  // Loop over the Q's and initialize each one.
  for (iterator QItr = this->begin(); QItr != this->end(); ++QItr) {
    (*QItr)->initialize(dataBase, boundaryBegin, boundaryEnd, time, dt, W);
  }
}

//------------------------------------------------------------------------------
// Method to calculate and return the viscous acceleration, work, and pressure,
// all in one step (efficiency and all).
//------------------------------------------------------------------------------
template<typename Dimension>
void
ArtificialViscosityList<Dimension>::
viscousEffects(typename Dimension::Vector& acceleration,
               typename Dimension::Scalar& work,
               typename Dimension::Scalar& pressure,
               const NodeIteratorBase<Dimension>& nodeI,
               const NodeIteratorBase<Dimension>& nodeJ,
               const typename Dimension::Vector& rij, 
               const typename Dimension::Vector& rijUnit,
               const typename Dimension::Vector& vi,
               const typename Dimension::Vector& vj,
               const typename Dimension::Vector& etai,
               const typename Dimension::Vector& etaj,
               const typename Dimension::Scalar ci,
               const typename Dimension::Scalar cj,
               const typename Dimension::Scalar Pi,
               const typename Dimension::Scalar Pj,
               const typename Dimension::Scalar rhoi,
               const typename Dimension::Scalar rhoj,
               const typename Dimension::Scalar hi,
               const typename Dimension::Scalar hj,
               const typename Dimension::Vector& gradW) const {
  cdebug << "ArtificialViscosityList::viscousEffects" << endl;

  REQUIRE(rhoi > 0.0);

  // Zero out the incoming data.
  acceleration.Zero();
  work = 0.0;
  pressure = 0.0;

  // Accumulate the acceleration, work, and pressure into a single values.
  Vector localAcceleration;
  Scalar localWork;
  Scalar localPressure;
  for (const_iterator QItr = this->begin(); QItr != this->end(); ++QItr) {
    cdebug << "Invoking viscousEffects on " << *QItr << endl;
    (*QItr)->viscousEffects(localAcceleration, localWork, localPressure,
                            nodeI, nodeJ, rij, rijUnit,
                            vi, vj, etai, etaj, ci, cj, Pi, Pj, rhoi, rhoj, 
                            hi, hj, gradW);
    acceleration += localAcceleration;
    work += localWork;
    pressure += localPressure;
  }
  cdebug << "Done." << endl;
}

//------------------------------------------------------------------------------
// Add an artificial viscosity to the list.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ArtificialViscosityList<Dimension>::
appendArtificialViscosity(ArtificialViscosity<Dimension>* QPtr) {
  // Add this Q to the list, if it's not already there.
  if (!haveArtificialViscosity(QPtr)) push_back(QPtr);
}

//------------------------------------------------------------------------------
// Delete an artificial viscosity from the list.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ArtificialViscosityList<Dimension>::
deleteArtificialViscosity(ArtificialViscosity<Dimension>* QPtr) {

  // Find an iterator pointing to this Q in the list.
  iterator QItr = find(this->begin(), this->end(), QPtr);
  
  // Remove the Q from the list, if we have it.
  if (QItr < this->end()) erase(QItr);
}

//------------------------------------------------------------------------------
// Test if the given artificial viscosity is listed in this list.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
ArtificialViscosityList<Dimension>::
haveArtificialViscosity(const ArtificialViscosity<Dimension>* QPtr) const {
  return find(this->begin(), this->end(), QPtr) != this->end();
}

//------------------------------------------------------------------------------
// Test if the ArtificialViscosity is valid, i.e., ready to use.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
ArtificialViscosityList<Dimension>::valid() const {
  bool ok = this->size() > 0;
  const_iterator QItr = this->begin();
  while (ok && QItr != this->end()) {
    ok = (*QItr)->valid();
    ++QItr;
  }
  return ok;
}
}
}

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
namespace Spheral {
namespace ArtificialViscositySpace {
template class ArtificialViscosityList< Dim<1> >;
template class ArtificialViscosityList< Dim<2> >;
template class ArtificialViscosityList< Dim<3> >;
}
}
