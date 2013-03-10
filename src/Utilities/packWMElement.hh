//---------------------------------Spheral++----------------------------------//
// packWMElement
// Methods to pack/unpack WildMagic objects into vector<char>.
//
// Created by J. Michael Owen, Mon Jan 25 20:36:56 PST 2010
//----------------------------------------------------------------------------//
#ifndef __Spheral_packWMElement__
#define __Spheral_packWMElement__

#include <vector>
#include "packElement.hh"
#include "Geometry/Dimension.hh"

#ifdef USE_MPI
extern "C" {
#include "mpi.h"
}
#endif

namespace Spheral {

//------------------------------------------------------------------------------
// Pack the WM Vectors.
//------------------------------------------------------------------------------
template<>
inline
void
packElement(const Dim<2>::WMVector& value,
            std::vector<char>& buffer) {
  packElement(value[0], buffer);
  packElement(value[1], buffer);
}

template<>
inline
void
packElement(const Dim<3>::WMVector& value,
            std::vector<char>& buffer) {
  packElement(value[0], buffer);
  packElement(value[1], buffer);
  packElement(value[2], buffer);
}

//------------------------------------------------------------------------------
// Pack the WM Boxes.
//------------------------------------------------------------------------------
template<>
inline
void
packElement(const Dim<2>::Box& value,
            std::vector<char>& buffer) {
  packElement(value.Center, buffer);
  packElement(value.Axis[0], buffer);
  packElement(value.Axis[1], buffer);
  packElement(value.Extent[0], buffer);
  packElement(value.Extent[1], buffer);
}

template<>
inline
void
packElement(const Dim<3>::Box& value,
            std::vector<char>& buffer) {
  packElement(value.Center, buffer);
  packElement(value.Axis[0], buffer);
  packElement(value.Axis[1], buffer);
  packElement(value.Axis[2], buffer);
  packElement(value.Extent[0], buffer);
  packElement(value.Extent[1], buffer);
  packElement(value.Extent[2], buffer);
}

//------------------------------------------------------------------------------
// Unpack the WM Vectors.
//------------------------------------------------------------------------------
template<>
inline
void
unpackElement(Dim<2>::WMVector& value,
              std::vector<char>::const_iterator& itr,
              const std::vector<char>::const_iterator& endPackedVector) {
  unpackElement(value[0], itr, endPackedVector);
  unpackElement(value[1], itr, endPackedVector);
}

template<>
inline
void
unpackElement(Dim<3>::WMVector& value,
              std::vector<char>::const_iterator& itr,
              const std::vector<char>::const_iterator& endPackedVector) {
  unpackElement(value[0], itr, endPackedVector);
  unpackElement(value[1], itr, endPackedVector);
  unpackElement(value[2], itr, endPackedVector);
}

//------------------------------------------------------------------------------
// Unpack the WM Boxes.
//------------------------------------------------------------------------------
template<>
inline
void
unpackElement(Dim<2>::Box& value,
              std::vector<char>::const_iterator& itr,
              const std::vector<char>::const_iterator& endPackedVector) {
  unpackElement(value.Center, itr, endPackedVector);
  unpackElement(value.Axis[0], itr, endPackedVector);
  unpackElement(value.Axis[1], itr, endPackedVector);
  unpackElement(value.Extent[0], itr, endPackedVector);
  unpackElement(value.Extent[1], itr, endPackedVector);
}

template<>
inline
void
unpackElement(Dim<3>::Box& value,
              std::vector<char>::const_iterator& itr,
              const std::vector<char>::const_iterator& endPackedVector) {
  unpackElement(value.Center, itr, endPackedVector);
  unpackElement(value.Axis[0], itr, endPackedVector);
  unpackElement(value.Axis[1], itr, endPackedVector);
  unpackElement(value.Axis[2], itr, endPackedVector);
  unpackElement(value.Extent[0], itr, endPackedVector);
  unpackElement(value.Extent[1], itr, endPackedVector);
  unpackElement(value.Extent[2], itr, endPackedVector);
}

}

#endif
