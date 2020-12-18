namespace Spheral {

//------------------------------------------------------------------------------
// Lookup the kernel for (rj/h, ri/h) = (etaj, etai)
//------------------------------------------------------------------------------
inline
double
SphericalTableKernel::operator()(const Dim<1>::Vector& etaj,
                                 const Dim<1>::Vector& etai) const {
  const auto ei = etai[0];
  const auto ej = etaj[0];
  CHECK(ei > 0.0);
  CHECK(ej > 0.0);
  return 2.0*M_PI/(ei*ej)*mInterp(Dim<2>::Vector(ej, ei));
}

//------------------------------------------------------------------------------
// kernelValue, gradValue, grad2Value -- pass through to the base TableKernel
//------------------------------------------------------------------------------
inline
double 
SphericalTableKernel::kernelValue(const double etaMagnitude, const double Hdet) const {
  return mKernel.kernelValue(etaMagnitude, Hdet);
}

inline
double 
SphericalTableKernel::gradValue(const double etaMagnitude, const double Hdet) const {
  return mKernel.gradValue(etaMagnitude, Hdet);
}

inline
double 
SphericalTableKernel::grad2Value(const double etaMagnitude, const double Hdet) const {
  return mKernel.grad2Value(etaMagnitude, Hdet);
}

}
