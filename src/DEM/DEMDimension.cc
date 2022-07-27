#include "DEMDimension.hh"

namespace Spheral{
    const Dim<1>::Scalar DEMDimension<Dim<1>>::zero = 0.0;
    const Dim<2>::Scalar DEMDimension<Dim<2>>::zero = 0.0;
    const Dim<3>::Vector DEMDimension<Dim<3>>::zero = Dim<3>::Vector(0.0, 0.0, 0.0);   
}