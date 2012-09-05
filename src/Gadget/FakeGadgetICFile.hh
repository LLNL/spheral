//---------------------------------Spheral++----------------------------------//
// FakeGadgetICFile -- An object that creates a bogus IC file for Gadget.
// This wouldn't be necessary if Gadget weren't L-A-M-E about how it reads its 
// initial conditions.
//
//! \author $Author: jeffjohnson $
//! \version $Revision: 434 $
//! \date $Date: 2002-10-11 00:36:55 -0700 (Fri, 11 Oct 2002) $
//----------------------------------------------------------------------------//
#ifndef FAKEGADGETICFILE_HH
#define FAKEGADGETICFILE_HH

#include "Geometry/Dimension.hh"
#include "DataBase/DataBase.hh"

#ifdef USE_GADGET

namespace Spheral {
namespace GadgetSpace {

using namespace Spheral::DataBaseSpace;

class FakeGadgetICFile {
public:

  //! Constructor.  This creates a fake IC file for use by Gadget.
  FakeGadgetICFile(const DataBase<Dim<3> >& db);

  //! Destructor.
  ~FakeGadgetICFile();

  //! Name of the IC file.
  const char* name() const;

  //! File format.
  int format() const;

private:

  // No default constructor, copy constructor or assignment operator.
  FakeGadgetICFile();
  FakeGadgetICFile(const FakeGadgetICFile&);
  FakeGadgetICFile& operator=(const FakeGadgetICFile&);

  // Filename of bogus IC file.
  char mFilename[100];

}; // end class GravityForce

} // end namespace GadgetSpace

} // end namespace Spheral

#endif // USE_GADGET

#endif
