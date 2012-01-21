//---------------------------------Spheral++----------------------------------//
// FluidNodeTraits -- A trait class to define traits specific to different
//   types of fluid nodes.
// 
// Created by JMO, Tue Oct 26 14:35:28 PDT 1999
//----------------------------------------------------------------------------//

#ifndef FLUIDNODETRAITS_HH
#define FLUIDNODETRAITS_HH

#include "Geometry/Dimension.hh"

template<typename Dimension> class SphNodeList;

template<typename Dimension, typename FluidNodeType> class FluidNodeTraits {};

template<>
class FluidNodeTraits<Dim<1>, SphNodeList<Dim<1> > > {
public:
  typedef Dim<1>::Scalar HType;
};

template<>
class FluidNodeTraits<Dim<2>, SphNodeList<Dim<2> > > {
public:
  typedef Dim<2>::Scalar HType;
};

template<>
class FluidNodeTraits<Dim<3>, SphNodeList<Dim<3> > > {
public:
  typedef Dim<3>::Scalar HType;
};

#endif
