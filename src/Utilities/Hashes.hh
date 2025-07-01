//---------------------------------Spheral++----------------------------------//
// hashes
//
// A collection of simple useful hashes.
//
// Created by JMO, Fri Jan  7 21:56:22 PST 2011
//----------------------------------------------------------------------------//
#ifndef __Spheral_hashes__
#define __Spheral_hashes__

namespace std {

template<>
struct hash<std::pair<size_t, size_t>> {
  size_t operator()(const std::pair<size_t, size_t>& p) const {
    return hash<size_t>()(p.first) ^ hash<size_t>()(p.second);
  }
};

} // namespace std

#endif

