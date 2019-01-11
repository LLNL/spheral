namespace Spheral {

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
inline
NinthOrderPolynomialFit::
NinthOrderPolynomialFit(const double C0,
                        const double C1, 
                        const double C2, 
                        const double C3, 
                        const double C4, 
                        const double C5, 
                        const double C6, 
                        const double C7, 
                        const double C8, 
                        const double C9):
  mC0(C0),
  mC1(C1),
  mC2(C2),
  mC3(C3),
  mC4(C4),
  mC5(C5),
  mC6(C6),
  mC7(C7),
  mC8(C8),
  mC9(C9) {
}

//------------------------------------------------------------------------------
// Copy constructor
//------------------------------------------------------------------------------
inline
NinthOrderPolynomialFit::
NinthOrderPolynomialFit(const NinthOrderPolynomialFit& rhs):
  mC0(rhs.mC0),
  mC1(rhs.mC1),
  mC2(rhs.mC2),
  mC3(rhs.mC3),
  mC4(rhs.mC4),
  mC5(rhs.mC5),
  mC6(rhs.mC6),
  mC7(rhs.mC7),
  mC8(rhs.mC8),
  mC9(rhs.mC9) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
inline
NinthOrderPolynomialFit::
~NinthOrderPolynomialFit() {
}

//------------------------------------------------------------------------------
// Assignment.
//------------------------------------------------------------------------------
inline
NinthOrderPolynomialFit&
NinthOrderPolynomialFit::
operator=(const NinthOrderPolynomialFit& rhs) {
  if (this != &rhs) {
    mC0 = rhs.mC0;
    mC1 = rhs.mC1;
    mC2 = rhs.mC2;
    mC3 = rhs.mC3;
    mC4 = rhs.mC4;
    mC5 = rhs.mC5;
    mC6 = rhs.mC6;
    mC7 = rhs.mC7;
    mC8 = rhs.mC8;
    mC9 = rhs.mC9;
  }
  return *this;
}

//------------------------------------------------------------------------------
// Compute the polynomial for the given argument.
//------------------------------------------------------------------------------
inline
double
NinthOrderPolynomialFit::
operator()(const double x) const {
  return (mC0 + 
          mC1*x +
          mC2*x*x +
          mC3*x*x*x +
          mC4*x*x*x*x +
          mC5*x*x*x*x*x +
          mC6*x*x*x*x*x*x +
          mC7*x*x*x*x*x*x*x +
          mC8*x*x*x*x*x*x*x*x +
          mC9*x*x*x*x*x*x*x*x*x);
}

}
