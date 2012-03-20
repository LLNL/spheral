//---------------------------------Spheral++----------------------------------//
//! \author $Author: jeffjohnson $
//! \version $Id: FakeGadgetICFile.cc 449 2002-10-26 06:57:28Z jeffjohnson $
//! \date $Date: 2002-10-25 23:57:28 -0700 (Fri, 25 Oct 2002) $
//----------------------------------------------------------------------------//

#ifdef USE_GADGET

#include "FakeGadgetICFile.hh"
#include "FieldList.hh"
#include <algorithm>

extern "C" {
#include "allvars.h" // Gadget's global variables.
#include "proto.h"   // Gadget's function prototypes.
}

//------------------------------------------------------------------------------
//                          Implementation notes
//------------------------------------------------------------------------------
// The format of the file was taken from the User's Guide to Serial Gadget, 
// in the section on Snapshot files.  There are a *LOT* of calls to this 
// macro SKIP, which reads an integer and sometimes expects that integer value 
// to be 256.  I have NO IDEA why this is, but it's easy enough to 
// accommodate.  All in the name of SPH, I suppose.
//
// Also, Gadget assumes (not surprisingly, in retrospect) that there is at
// least one particle in the system, and does not test for zero-length 
// particle lists.  So we have to tell it to use at least one particle in our 
// IC file.  Beware nasty C codes!
//------------------------------------------------------------------------------

namespace Spheral {
namespace GadgetSpace {

using namespace Spheral::FieldSpace;
  
//------------------------------------------------------------------------------
FakeGadgetICFile::FakeGadgetICFile(const DataBase<Dim<3> >& db)
{
  // Some type definitions for convenience.
  typedef Dim<3>::Scalar Scalar;
  typedef Dim<3>::Vector Vector;

  // Count the number of particles in Spheral.
  size_t numberOfSpheralParticles = 0;
  for (DataBase<Dim<3> >::ConstNodeListIterator 
      nodeListIter = db.nodeListBegin();
      nodeListIter != db.nodeListEnd();
      ++nodeListIter) {
    numberOfSpheralParticles += (*nodeListIter)->numInternalNodes();
  } // end for

  // Create a temporary file and get a file descriptor for it.
  strcpy(mFilename, "spheral-gadget-XXXXXX");
  int bogusFileDescriptor = mkstemp(mFilename);
  FILE* bogusFile = fdopen(bogusFileDescriptor, "w");

  //------------------------------------------------------------------
  //                        Generate the header:
  //------------------------------------------------------------------

  // There's a really annoying macro in read_ic.c of Gadget that reads a 
  // single integer from the file.  This is called SKIP and is used at 
  // various points in the file to separate sections, as far as I can tell.
  // For no apparent reason, sometimes Gadget expects to read the value 
  // 256.  So I'm using my own SKIP macro to write out these "spacers."
  // This macro will be #undef'd at the end of this ordeal, lest you wish
  // to come after me with a large knife.
  int twoFiftySix = 256;
  #define SKIP std::fwrite(reinterpret_cast<void*>(&twoFiftySix), sizeof(int), 1, bogusFile)

  // It's convenient to keep some arrays of zeros around.
  double floatZeros[20];
  std::fill(floatZeros, floatZeros + 20, 0.0);
  int intZeros[20];
  std::fill(intZeros, intZeros + 20, 0);

  // Here we go!
  SKIP;
  
  // We are going to use Gadget halo particles to replicate Spheral's 
  // particles.  Why halo particles?  Because gas particles receive 
  // hydrodynamic acceleration.
  int numbersOfParticles[6];
  numbersOfParticles[0] = 0;
  numbersOfParticles[1] = numberOfSpheralParticles;
  std::fill(numbersOfParticles + 2, numbersOfParticles + 6, 0);
  std::fwrite(reinterpret_cast<void*>(numbersOfParticles), sizeof(int), 6, bogusFile);

  // Particle masses by type.  These are set to zero so that the particle 
  // masses will be read directly from Spheral's database.
  double particleTypeMasses[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  std::fwrite(reinterpret_cast<void*>(particleTypeMasses), sizeof(double), 6, bogusFile);

  // Zeros for the time expansion factor and the redshift, respectively.
  std::fwrite(reinterpret_cast<void*>(floatZeros), sizeof(double), 2, bogusFile);

  // Zeros for the star formation and feedback options, which are 
  // unavailable in the public version of Gadget.
  std::fwrite(reinterpret_cast<void*>(intZeros), sizeof(int), 2, bogusFile);

  // In the serial version, write out the numbers of particle types again.
  std::fwrite(reinterpret_cast<void*>(numbersOfParticles), sizeof(int), 6, bogusFile);

  // Disable cooling.
  std::fwrite(reinterpret_cast<void*>(intZeros), sizeof(int), 1, bogusFile);

  // The number of files containing this snapshot.
  int one = 1;
  std::fwrite(reinterpret_cast<void*>(&one), sizeof(int), 1, bogusFile);

  // Zero out the cosmological parameters.
  std::fwrite(reinterpret_cast<void*>(floatZeros), sizeof(double), 4, bogusFile);

  // Now write out the rest of the header, which is the remainder of 256 
  // bytes.  By my calculation, 
  // remainder = 256 - 6*4 - 6*8 - 2*8 - 2*4 - 6*4 - 2*4 - 4*8 = 96.
  // These can all be zeros.
  char padding[96];
  std::fill(padding, padding + 96, 0);
  std::fwrite(reinterpret_cast<void*>(padding), sizeof(char), 96, bogusFile);

  // Gadget wants to skip again.
  SKIP;

  // Now we need to go to the database and get information about our particles.

  // Allocate memory for the particle information transfer.
  float* particleFloats = new float[3*numberOfSpheralParticles];
  int* particleInts = new int[numberOfSpheralParticles];

  // Position 
  const FieldList<Dim<3>, Vector>& position = 
    db.globalPosition();
  size_t i = 0;
  for (DataBase<Dim<3> >::IDIterator 
       nodeIter = db.internalNodeBegin();
       nodeIter != db.internalNodeEnd();
       ++nodeIter) {
    particleFloats[3*i] = static_cast<float>(position(nodeIter)(0));
    particleFloats[3*i+1] = static_cast<float>(position(nodeIter)(1));
    particleFloats[3*i+2] = static_cast<float>(position(nodeIter)(2));
    ++i;
  } // end for
  SKIP;
  std::fwrite(reinterpret_cast<void*>(particleFloats), sizeof(float), 
              3*numberOfSpheralParticles, bogusFile);
  SKIP;

  
  // Velocity 
  const FieldList<Dim<3>, Vector>& velocity = 
    db.globalVelocity();
  i = 0;
  for (DataBase<Dim<3> >::IDIterator 
       nodeIter = db.internalNodeBegin();
       nodeIter != db.internalNodeEnd();
       ++nodeIter) {
    particleFloats[3*i] = static_cast<float>(velocity(nodeIter)(0));
    particleFloats[3*i+1] = static_cast<float>(velocity(nodeIter)(1));
    particleFloats[3*i+2] = static_cast<float>(velocity(nodeIter)(2));
    ++i;
  } // end for
  SKIP;
  std::fwrite(reinterpret_cast<void*>(particleFloats), sizeof(float), 
              3*numberOfSpheralParticles, bogusFile);
  SKIP;
 
  // Particle ID.
  for (i = 0; i < numberOfSpheralParticles; ++i)
  {
    particleInts[i] = i+1;
  } // end for
  SKIP;
  std::fwrite(reinterpret_cast<void*>(particleInts), sizeof(int), 
              numberOfSpheralParticles, bogusFile);
  SKIP;

  // Particle masses.
  const FieldList<Dim<3>, Scalar>& mass = 
    db.globalMass();
  i = 0;
  for (DataBase<Dim<3> >::IDIterator 
       nodeIter = db.internalNodeBegin();
       nodeIter != db.internalNodeEnd();
       ++nodeIter) {
      particleFloats[i] = static_cast<float>(mass(nodeIter));
      ++i;
    } // end for
    SKIP;
    std::fwrite(reinterpret_cast<void*>(particleFloats), sizeof(float), 
        numberOfSpheralParticles, bogusFile);
    SKIP;

    // Particle specific energy (none for Gadget's purposes).
    SKIP;
    std::fill(particleFloats, particleFloats + numberOfSpheralParticles, 0.0);
    std::fwrite(reinterpret_cast<void*>(particleFloats), sizeof(float), 
        numberOfSpheralParticles, bogusFile);

    SKIP;
/*
    // Particle density for comoving SPH or something (0).
    std::fwrite(reinterpret_cast<void*>(particleFloats), sizeof(float), 
        numberOfSpheralParticles, bogusFile);
    SKIP;
*/

    // Clean up the mess we've made.  (I hope we're happy with ourselves!!!)
    delete [] particleFloats;
    delete [] particleInts;
#undef SKIP

    // Now we're done.
    fclose(bogusFile);

}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
FakeGadgetICFile::
~FakeGadgetICFile() {
  // Get rid of this lame-ass file.
  int status = remove(mFilename);
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
const char*
FakeGadgetICFile::
name() const {
  // Get rid of this lame-ass file.
  return const_cast<const char*>(mFilename);
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
int
FakeGadgetICFile::
format() const {
  return 1;
}
//------------------------------------------------------------------------------

} // end namespace GadgetSpace
} // end namespace Spheral

#endif // USE_GADGET
