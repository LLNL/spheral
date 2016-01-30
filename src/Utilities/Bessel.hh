#ifndef SPHERAL_BESSEL_HH
#define SPHERAL_BESSEL_HH

namespace Spheral
{

//! \namespace Bessel
//! The functions contained in this namespace are the Bessel functions of 
//! integer and real order.
namespace Bessel
{

//! Compute J0(x), the Bessel function of the first kind of order 0.
//! \param x The argument to the Bessel function.
double J0(double x);

//! Compute J1(x), the Bessel function of the first kind of order 1.
//! \param x The argument to the Bessel function.
double J1(double x);

//! Compute Y0(x), the Bessel function of the second kind of order 0.
double Y0(double x);

//! Compute Y1(x), the Bessel function of the second kind of order 1.
double Y1(double x);

//! Compute the nth zero of J0.
//! \param n The zero of J0 to be computed.
double zeroOfJ0(unsigned int n);

//! Compute the nth zero of J1.
//! \param n The zero of J1 to be computed.
double zeroOfJ1(unsigned int n);

}
}

#endif
