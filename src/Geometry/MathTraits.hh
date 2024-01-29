//---------------------------------Spheral++----------------------------------//
// MathTraits -- trait class for dimensional type information.
//
// Created by JMO, Fri Apr 16 10:13:52 PDT 1999
// -- JMO, Sun Apr 30 13:17:31 PDT 2000
//    replacing old defunct class with new trait class to define div and
//    grad types.
//----------------------------------------------------------------------------//

#ifndef MathTraits_HH
#define MathTraits_HH

#include "Dimension.hh"

namespace Spheral {

//------------------------------------------------------------------------------
template<typename Dimension, typename DataType> class MathTraits {};

//------------------------------------------------------------------------------
template<>
class MathTraits<Dim<1>, Dim<1>::Scalar> {
public:
  typedef Dim<1>::Vector GradientType;
  typedef Dim<1>::Scalar DivergenceType;
  typedef Dim<1>::Tensor HessianType;
};

template<>
class MathTraits<Dim<1>, Dim<1>::Vector> {
public:
  typedef Dim<1>::Tensor GradientType;
  typedef Dim<1>::Scalar DivergenceType;
  typedef Dim<1>::ThirdRankTensor HessianType;
};

template<>
class MathTraits<Dim<1>, Dim<1>::Tensor> {
public:
  typedef Dim<1>::Tensor GradientType;
  typedef Dim<1>::Vector DivergenceType;
};

template<>
class MathTraits<Dim<1>, Dim<1>::SymTensor> {
public:
  typedef Dim<1>::Tensor GradientType;
  typedef Dim<1>::Vector DivergenceType;
};

//------------------------------------------------------------------------------
template<>
class MathTraits<Dim<2>, Dim<2>::Scalar> {
public:
  typedef Dim<2>::Vector GradientType;
  typedef Dim<2>::Scalar DivergenceType;
  typedef Dim<2>::Tensor HessianType;
};

template<>
class MathTraits<Dim<2>, Dim<2>::Vector> {
public:
  typedef Dim<2>::Tensor GradientType;
  typedef Dim<2>::Scalar DivergenceType;
  typedef Dim<2>::ThirdRankTensor HessianType;
};

template<>
class MathTraits<Dim<2>, Dim<2>::Tensor> {
public:
  typedef Dim<2>::Tensor GradientType;
  typedef Dim<2>::Vector DivergenceType;
};

template<>
class MathTraits<Dim<2>, Dim<2>::SymTensor> {
public:
  typedef Dim<2>::Tensor GradientType;
  typedef Dim<2>::Vector DivergenceType;
};

//------------------------------------------------------------------------------
template<>
class MathTraits<Dim<3>, Dim<3>::Scalar> {
public:
  typedef Dim<3>::Vector GradientType;
  typedef Dim<3>::Scalar DivergenceType;
  typedef Dim<3>::Tensor HessianType;
};

template<>
class MathTraits<Dim<3>, Dim<3>::Vector> {
public:
  typedef Dim<3>::Tensor GradientType;
  typedef Dim<3>::Scalar DivergenceType;
  typedef Dim<3>::ThirdRankTensor HessianType;
};

template<>
class MathTraits<Dim<3>, Dim<3>::Tensor> {
public:
  typedef Dim<3>::Tensor GradientType;
  typedef Dim<3>::Vector DivergenceType;
};

template<>
class MathTraits<Dim<3>, Dim<3>::SymTensor> {
public:
  typedef Dim<3>::Tensor GradientType;
  typedef Dim<3>::Vector DivergenceType;
};
}

#endif
