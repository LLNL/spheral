//---------------------------------Spheral++----------------------------------//
// GeomFifthRankTensor -- Geometric Tensor (rank 5) Class.
//
// This is a very simple, limited functionality rank 5 tensor.  Assumes
// the thing is dimensioned nDim x nDim x nDim x nDim.
//
// Created by JMO, Tue Oct 13 15:24:56 PDT 2015
//----------------------------------------------------------------------------//
#ifndef __Spheral_GeomFifthRankTensor_hh__
#define __Spheral_GeomFifthRankTensor_hh__

#include <iostream>

#include "Geometry/RankNTensor.hh"
#include "Geometry/GeomFifthRankTensor_fwd.hh"

namespace Spheral {

template<int nDim>
class GeomFifthRankTensor : public RankNTensor<nDim, 5, GeomFifthRankTensor<nDim> > {
  using BaseType = RankNTensor<nDim, 5, GeomFifthRankTensor<nDim>>;

public:
  //--------------------------- Public Interface ---------------------------//
  using size_type = typename BaseType::size_type;
  static constexpr size_type numElements = BaseType::numElements;

  // Useful static member data.
  static const GeomFifthRankTensor zero;

  // Constructors.
  SPHERAL_HOST_DEVICE GeomFifthRankTensor();
  SPHERAL_HOST_DEVICE explicit GeomFifthRankTensor(const double val);
  SPHERAL_HOST_DEVICE GeomFifthRankTensor(const GeomFifthRankTensor& rhs);

  // Destructor.
  SPHERAL_HOST_DEVICE ~GeomFifthRankTensor();

  // Assignment.
  SPHERAL_HOST_DEVICE GeomFifthRankTensor& operator=(const GeomFifthRankTensor& rhs);
  SPHERAL_HOST_DEVICE GeomFifthRankTensor& operator=(const double rhs);

  // Access the elements by indicies.
  SPHERAL_HOST_DEVICE double operator()(const size_type i, const size_type j, const size_type k, const size_type m, const size_type n) const;
  SPHERAL_HOST_DEVICE double& operator()(const size_type i, const size_type j, const size_type k, const size_type m, const size_type n);

private:
  //--------------------------- Private Interface ---------------------------//
  using RankNTensor<nDim, 5, GeomFifthRankTensor>::mElements;
};

template<int nDims> const GeomFifthRankTensor<nDims> GeomFifthRankTensor<nDims>::zero = GeomFifthRankTensor<nDims>(0.0);

}

#include "GeomFifthRankTensorInline.hh"

#endif
