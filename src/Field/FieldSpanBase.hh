//---------------------------------Spheral++----------------------------------//
// FieldSpanBase -- An abstract base class to provide a generic handle on
// FieldSpans.
//
// See FieldSpan for more detailed information and the relation to Fields
//
// Created by JMO, Mon Apr 28 14:25:28 PDT 2025
//----------------------------------------------------------------------------//
#ifndef __Spheral_FieldSpanBase__
#define __Spheral_FieldSpanBase__

namespace Spheral {

template<typename Dimension>
class FieldSpanBase {

public:
  //--------------------------- Public Interface ---------------------------//
  // Constructors, destructor
  FieldSpanBase() = default;
  FieldSpanBase(const FieldSpanBase& rhs) = default;
  virtual ~FieldSpanBase() = default;

  // Assignment operator.
  virtual FieldSpanBase& operator=(const FieldSpanBase& rhs) = 0;

  // Require descendent fields be able to test equivalence.
  virtual bool operator==(const FieldSpanBase& rhs) const = 0;
  bool operator!=(const FieldSpanBase& rhs) const { return !(*this == rhs); }

  // Methods every FieldSpan must provide.
  virtual size_t size() const = 0;
  virtual void Zero() = 0;
};

}

#endif
