//---------------------------------Spheral++----------------------------------//
// PolynomialFit -- A ninth order polynomial fitting class for use with the
// Steinberg-Guninan strenth model.
//
// Created by JMO, Thu Sep 9 17:02:38 2004
//----------------------------------------------------------------------------//
#ifndef __Spheral_PolynomialFit_hh__
#define __Spheral_PolynomialFit_hh__

namespace Spheral {

class NinthOrderPolynomialFit {
public:
  //--------------------------- Public Interface ---------------------------//
  // Constructors, destructor.
  NinthOrderPolynomialFit(const double C0, 
                          const double C1, 
                          const double C2, 
                          const double C3, 
                          const double C4, 
                          const double C5, 
                          const double C6, 
                          const double C7, 
                          const double C8, 
                          const double C9);
  NinthOrderPolynomialFit(const NinthOrderPolynomialFit& rhs);
  ~NinthOrderPolynomialFit();
                
  // Assignment.
  NinthOrderPolynomialFit& operator=(const NinthOrderPolynomialFit& rhs);

  // Compute the polynomial for the given argument.
  double operator()(const double x) const;

private:
  //--------------------------- Private Interface ---------------------------//
  double mC0;
  double mC1;
  double mC2;
  double mC3;
  double mC4;
  double mC5;
  double mC6;
  double mC7;
  double mC8;
  double mC9;

  // No default constructor.
  NinthOrderPolynomialFit();
};

}

#include "PolynomialFitInline.hh"

#endif

