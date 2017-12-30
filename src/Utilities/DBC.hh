#ifndef DBC_HH
#define DBC_HH
//---------------------------------------------------------------------------
//
// DBC.hh -- Design by Contract tools.
//
//---------------------------------------------------------------------------
#include <iostream>
#include <sstream>
#include <string>
#include <exception>
#include <cmath>
#include "Process.hh"

#ifndef DBC_FUNCTIONS_HH
#define DBC_FUNCTIONS_HH

namespace Spheral
{
namespace dbc
{

// Lock for assertions
bool assertionLock (void);
void assertionUnLock(void);

// Exception thrown if VERIFY fails
class VERIFYError: public std::exception {
   public:

   // Constructor
   VERIFYError (const char* text):
      mReason(text) {}

   VERIFYError (const std::string& text):
      mReason(text) {}

   // Destructor
   ~VERIFYError() throw() {}

   // Description of the exception.
   const char* what() const throw() {
      return mReason.c_str();
   } // end what

   private:

   // Explanation of the failure.
   std::string mReason;
};


// nearlyEqual returns true if x and y approximately equal
template<typename T, typename U>
inline bool nearlyEqual(const T& x, 
                 const U& y, 
                 double relativeTolerance = 1.0e-5,
                 double absoluteTolerance = 1.0e-12)
{
   return (std::abs(x-y) <=
           (std::abs(y) * relativeTolerance + absoluteTolerance)
          );
}
} // namespace dbc
}

#endif // DBC_FUNCTIONS_HH

//----------------------------------------------------------------------------
//                         Clear any existing DBC compile flags.
//----------------------------------------------------------------------------

#ifdef DBC_USE_REQUIRE
#undef DBC_USE_REQUIRE
#endif

#ifdef DBC_USE_ENSURE
#undef DBC_USE_ENSURE
#endif

#ifdef DBC_USE_SCOPES
#undef DBC_USE_SCOPES
#endif

//----------------------------------------------------------------------------
//                         Work out the compilation modes.
//----------------------------------------------------------------------------

#ifdef DBC_COMPILE_ALL
#define DBC_USE_REQUIRE
#define DBC_USE_ENSURE
#define DBC_USE_SCOPES
#endif

#ifdef DBC_COMPILE_PRE
#define DBC_USE_REQUIRE
#define DBC_USE_SCOPES
#endif

//----------------------------------------------------------------------------
//                            REQUIRE -- Preconditions
//----------------------------------------------------------------------------

#ifdef DBC_USE_REQUIRE
#define DBC_ASSERTION(x, msg, kind)                     \
   if (::Spheral::dbc::assertionLock()) {               \
      if (!(x)) {                                       \
      std::stringstream s_SS;                           \
      s_SS << kind << ": " << msg << std::endl;         \
      s_SS << "...at line " << __LINE__ <<              \
         " of file " << __FILE__ << "." << std::ends;   \
      ::Spheral::Process::haltAll(s_SS.str().c_str());  \
   }                                                    \
   ::Spheral::dbc::assertionUnLock();                   \
}
#define REQUIRE2(x, msg) DBC_ASSERTION(x, msg, "Precondition violated")
#define ASSERT2(x, msg) DBC_ASSERTION(x, msg, "Assertion violated")
#else
#define ASSERT2(x, msg)
#define REQUIRE2(x, msg)
#endif

//----------------------------------------------------------------------------
//                            ENSURE -- Postconditions and Invariants
//----------------------------------------------------------------------------

//----------------------------------------------------------------------------
#ifdef DBC_USE_ENSURE
#define ENSURE2(x, msg) DBC_ASSERTION(x, msg,"Postcondition violated")
#define INVARIANT2(x, msg) DBC_ASSERTION(x, msg, "Invariant violated")
#else
#define ENSURE2(x, msg)
#define INVARIANT2(x, msg)
#endif
//----------------------------------------------------------------------------
//                               Contract scopes
// Use these to declare variables that are only used by contracts.
//----------------------------------------------------------------------------
#ifdef DBC_USE_SCOPES
#define BEGIN_CONTRACT_SCOPE {
#define END_CONTRACT_SCOPE }
#else
#define BEGIN_CONTRACT_SCOPE if (false) { while(false)
#define END_CONTRACT_SCOPE } while(false);
#endif

//----------- Define one-argument forms
#ifdef ASSERT
#undef ASSERT
#endif
#define ASSERT(x) ASSERT2(x, #x)
#define REQUIRE(x) REQUIRE2(x, #x)
#define ENSURE(x) ENSURE2(x, #x)
#define INVARIANT(x) INVARIANT2(x, #x)

//----------- Define alternate forms
#define CHECK(x) ASSERT(x)
#define CHECK2(x, msg) ASSERT2(x, msg)

#define VERIFY2(x, msg) \
   if (!(x)) { \
      std::stringstream s; \
      s << "Verification failed: " << msg << std::endl; \
      s << "...at line " << __LINE__ << \
         " of file " << __FILE__ << "." << std::ends;\
      ::Spheral::dbc::VERIFYError reason(s.str());\
      throw reason;\
   }
#define VERIFY(x) VERIFY2(x, #x)

// //----------------------------------------------------------------------------
// // Make lower case versions of all the contracts.
// //----------------------------------------------------------------------------
// #define require    REQUIRE
// #define ensure     ENSURE
// #define check      CHECK
// #define invariant  INVARIANT

// #define require2    REQUIRE2
// #define ensure2     ENSURE2
// #define check2      CHECK2
// #define invariant2  INVARIANT2

#endif // DBC_HH
