//---------------------------------Spheral++----------------------------------//
// GeomThirdRankTensor -- Geometric Tensor (rank 3) Class.
//
// This is a very simple, limited functionality rank 3 tensor.  Assumes
// the thing is dimensioned nDim x nDim x nDim.
//
// Created by JMO, Thu Jul 17 17:05:58 PDT 2008
//----------------------------------------------------------------------------//
#ifndef __Spheral_GeomThirdRankTensor_hh__
#define __Spheral_GeomThirdRankTensor_hh__

#include <iostream>

#include "Geometry/RankNTensor.hh"
#include "Geometry/GeomThirdRankTensor_fwd.hh"

namespace Spheral {

template<int nDim>
class GeomThirdRankTensor : public RankNTensor<nDim, 3, GeomThirdRankTensor<nDim> > {
  using BaseType = RankNTensor<nDim, 3, GeomThirdRankTensor<nDim>>;

public:
  //--------------------------- Public Interface ---------------------------//
  using size_type = typename BaseType::size_type;
  static constexpr size_type numElements = BaseType::numElements;

  // Useful static member data.
  static const GeomThirdRankTensor zero;

  // Constructors.
  GeomThirdRankTensor();
  explicit GeomThirdRankTensor(const double val);
  GeomThirdRankTensor(const GeomThirdRankTensor& rhs);

  // Destructor.
  ~GeomThirdRankTensor();

  // Assignment.
  GeomThirdRankTensor& operator=(const GeomThirdRankTensor& rhs);
  GeomThirdRankTensor& operator=(const double rhs);

  // Access the elements by indicies.
  double operator()(const size_type i, const size_type j, const size_type k) const;
  double& operator()(const size_type i, const size_type j, const size_type k);

private:
  //--------------------------- Private Interface ---------------------------//
  using RankNTensor<nDim, 3, GeomThirdRankTensor>::mElements;
};


template<int nDims> const GeomThirdRankTensor<nDims> GeomThirdRankTensor<nDims>::zero = GeomThirdRankTensor<nDims>(0.0);

}

#include "GeomThirdRankTensorInline.hh"

#endif
