//---------------------------------Spheral++----------------------------------//
// FieldSpanListBase -- A base class to provide a generic handle on
//                      FieldSpanLists.
//
// Created by JMO, Thu May  1 15:21:54 PDT 2025
//----------------------------------------------------------------------------//
#ifndef __Spheral_FieldSpanListBase__
#define __Spheral_FieldSpanListBase__

namespace Spheral {

template<typename Dimension>
class FieldSpanListBase {

public:
  //--------------------------- Public Interface ---------------------------//
  // Constructors, destructor
  FieldSpanListBase() = default;
  FieldSpanListBase(FieldSpanListBase& rhs) = default;
  FieldSpanListBase& operator=(const FieldSpanListBase& rhs) = default;
  virtual ~FieldSpanListBase() = default;
};

}

#endif
