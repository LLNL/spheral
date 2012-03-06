#ifndef __Spheral_SuffixTraits__
#define __Spheral_SuffixTraits__

#include <string>
#include "Geometry/Dimension.hh"

// Trait class to enforce the Spheral standard for dimension suffixes in 
// Python.
namespace Spheral {
  template<typename Dimension> struct SuffixTraits { static std::string suffix() {} };
  template<> struct SuffixTraits<Dim<1> > { static std::string suffix() { return "1d"; } };
  template<> struct SuffixTraits<Dim<2> > { static std::string suffix() { return "2d"; } };
  template<> struct SuffixTraits<Dim<3> > { static std::string suffix() { return "3d"; } };
}

#else

namespace Spheral {
  template<typename Dimension> struct SuffixTraits;
}

#endif
