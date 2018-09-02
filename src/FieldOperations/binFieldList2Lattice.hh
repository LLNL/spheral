//---------------------------------Spheral++----------------------------------//
// binFieldList2Lattice
//
// Bin the values of a FieldList to a lattice.
// The results are returned as a vector<Value>, of size nsample[0]*nsample[1]*....
// Note, in parallel this does the global reduction for you, so you get back
// the global result.
//
// Created by JMO, Wed Feb  3 11:02:20 PST 2010
//----------------------------------------------------------------------------//
#ifndef __Spheral_binField2Lattice__
#define __Spheral_binField2Lattice__

#include <vector>

namespace Spheral {

template<typename Dimension> class TableKernel;

template<typename Dimension, typename DataType> class FieldList;

//------------------------------------------------------------------------------
// Do a straightforward binning.
//------------------------------------------------------------------------------
template<typename Dimension, typename Value>
std::vector<Value>
binFieldList2Lattice(const FieldList<Dimension, Value>& fieldList,
                     const typename Dimension::Vector& xmin,
                     const typename Dimension::Vector& xmax,
                     const std::vector<unsigned>& nsample);

//------------------------------------------------------------------------------
// Bin with kernel smoothing applied to the FieldList.
//------------------------------------------------------------------------------
template<typename Dimension, typename Value>
std::vector<Value>
binFieldList2Lattice(const FieldList<Dimension, Value>& fieldList,
                     const TableKernel<Dimension>& W,
                     const typename Dimension::Vector& xmin,
                     const typename Dimension::Vector& xmax,
                     const std::vector<unsigned>& nsample);

//------------------------------------------------------------------------------
// Helper method to translate the given position to an index into a 1-D array.
//------------------------------------------------------------------------------
inline
size_t
latticeIndexElement(const double xi,
                    const double xmin,
                    const double xmax,
                    const unsigned nsample) {
  return std::min(nsample - 1U, unsigned(std::max(0.0, std::min(1.0, (xi - xmin)/(xmax - xmin)))*nsample));
}

inline
size_t
latticeIndex(const Dim<1>::Vector& xi,
             const Dim<1>::Vector& xmin,
             const Dim<1>::Vector& xmax,
             const std::vector<unsigned>& nsample) {
  REQUIRE(nsample.size() == 1);
  return latticeIndexElement(xi.x(), xmin.x(), xmax.x(), nsample[0]);
}

inline
size_t
latticeIndex(const Dim<2>::Vector& xi,
             const Dim<2>::Vector& xmin,
             const Dim<2>::Vector& xmax,
             const std::vector<unsigned>& nsample) {
  REQUIRE(nsample.size() == 2);
  return (latticeIndexElement(xi.x(), xmin.x(), xmax.x(), nsample[0]) +
          latticeIndexElement(xi.y(), xmin.y(), xmax.y(), nsample[1])*nsample[0]);
}

inline
size_t
latticeIndex(const Dim<3>::Vector& xi,
             const Dim<3>::Vector& xmin,
             const Dim<3>::Vector& xmax,
             const std::vector<unsigned>& nsample) {
  REQUIRE(nsample.size() == 3);
  return (latticeIndexElement(xi.x(), xmin.x(), xmax.x(), nsample[0]) +
          latticeIndexElement(xi.y(), xmin.y(), xmax.y(), nsample[1])*nsample[0] +
          latticeIndexElement(xi.z(), xmin.z(), xmax.z(), nsample[2])*nsample[0]*nsample[1]);
}

}

#endif
