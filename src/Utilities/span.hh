//---------------------------------Spheral++----------------------------------//
// Wrapper for std::span
//
// If we're using a C++20 compiler we can use std::span, otherwise use
// boost::span
//----------------------------------------------------------------------------//
#ifndef __Spheral_Span__
#define __Spheral_Span__

// Check for C++20 support
// #ifdef 1 __cpp_lib_span
  #define SPHERAL_SPAN_HEADER <span>
  #define SPHERAL_SPAN_TYPE std::span
// #else
//   #define SPHERAL_SPAN_HEADER <boost/core/span.hpp>
//   #define SPHERAL_SPAN_TYPE boost::span
// #endif

// Include the appropriate header
#include SPHERAL_SPAN_HEADER

#endif
