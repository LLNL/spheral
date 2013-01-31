//---------------------------------------------------------------------------
//
// Process.cc
//
//---------------------------------------------------------------------------

#include <cstdlib>
#include <iostream>

#include "Process.hh"
#include "Distributed/Communicator.hh"

#ifdef USE_MPI
#include <mpi.h>

// Static member data initialization.
int Spheral::Process::sRank = -1;
int Spheral::Process::sTotalProcs = -1;
#endif

namespace Spheral {
//----------------------------------------------------------------------------
int 
Process::
getRank()
{
#ifdef USE_MPI
   if (sRank == -1)
   {
      int isInitialized;
      MPI_Initialized(&isInitialized);
      if (isInitialized)
      {
         MPI_Comm_rank(Communicator::communicator(), &sRank);
      } // end if
      else
      {
         sRank = 0;
      } // end else
   } // end if
   return sRank;
#else
   return 0;
#endif
} // end getRank
//----------------------------------------------------------------------------


//----------------------------------------------------------------------------
int
Process::
getTotalNumberOfProcesses()
{
#ifdef USE_MPI
   if (sTotalProcs == -1)
   {
      int isInitialized;
      MPI_Initialized(&isInitialized);
      if (isInitialized)
      {
         MPI_Comm_size(Communicator::communicator(), &sTotalProcs);
      } // end if
      else
      {
         sTotalProcs = 1;
      } // end else
   } // end if
   return sTotalProcs;
#else
   return 1;
#endif
} // end getTotalNumberOfProcesses

void Process::haltAll(const char* msg) {
   std::cout << msg << std::endl;
   std::cout.flush();
   std::cerr.flush();  
#ifdef USE_MPI
   int isInitialized;
   MPI_Initialized(&isInitialized);
   if (isInitialized)
   {
      MPI_Abort(Communicator::communicator(), 1);
   } 
#endif
   abort();
} // end Process::haltAll
}
//----------------------------------------------------------------------------

