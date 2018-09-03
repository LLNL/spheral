#include "Utilities/Bessel.hh"
#include "Utilities/DBC.hh"
#include <algorithm>

// Use the GNU scientific library to evaluate the Bessel functions.
#include "gsl/gsl_sf_bessel.h"
#include "gsl/gsl_machine.h" // For GSL_DBL_EPSILON, etc.


namespace Spheral {
namespace Bessel {




//----------------------------------------------------------------------------
double 
J0(double x)
{
   return gsl_sf_bessel_J0(x);
} // end J0
//----------------------------------------------------------------------------


//----------------------------------------------------------------------------
double
J1(double x)
{
   return gsl_sf_bessel_J1(x); 
} // end J1
//----------------------------------------------------------------------------


//----------------------------------------------------------------------------
double
Y0(double x)
{
   REQUIRE2(x > 0.0, "Y0(x) is not defined for non-positive x!");
   // Fix the bounds to work with GSL if necessary.
   x = min(x, (0.99 / GSL_DBL_EPSILON)); // top out
   return gsl_sf_bessel_Y0(x);
} // end Y0
//----------------------------------------------------------------------------


//----------------------------------------------------------------------------
double
Y1(double x)
{
   REQUIRE2(x > 0.0, "Y1(x) is not defined for non-positive x!");
   // Fix the bounds to work with GSL if necessary.
   x = min(x, (0.99 / GSL_DBL_EPSILON)); // top out
   x = max(x, 2.01 * GSL_SQRT_DBL_EPSILON); // bottom out
   return gsl_sf_bessel_Y1(x);
} // end Y1
//----------------------------------------------------------------------------


//----------------------------------------------------------------------------
double
zeroOfJ0(unsigned int n)
{
   return gsl_sf_bessel_zero_J0(n);
} // end zeroOfJ0
//----------------------------------------------------------------------------


//----------------------------------------------------------------------------
double
zeroOfJ1(unsigned int n)
{
   return gsl_sf_bessel_zero_J1(n);
} // end zeroOfJ1
//----------------------------------------------------------------------------


}
}
