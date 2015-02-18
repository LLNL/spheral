//---------------------------------Spheral++----------------------------------//
// FieldDataTypeTraits, specializations of the DataTypeTraits for the Field
// class.
//
// Created by J. Michael Owen, Mon Nov 21 21:12:11 PST 2011
//----------------------------------------------------------------------------//
#ifndef FieldDataTypeTraits_HH
#define FieldDataTypeTraits_HH

#include "DataTypeTraits.hh"
#include "Field/Field.hh"

namespace Spheral {

//------------------------------------------------------------------------------
template<typename Dimension, typename Value>
struct DataTypeTraits<FieldSpace::Field<Dimension, Value> > {
  typedef Value ElementType;
  static bool fixedSize() { return false; }
  static int numElements(const FieldSpace::Field<Dimension, Value>& x) { return x.size(); }
  static FieldSpace::Field<Dimension, Value> zero() { return FieldSpace::Field<Dimension, Value>("Unnamed field"); }
};

//------------------------------------------------------------------------------
template<typename Dimension, typename Value>
struct DataTypeTraits<FieldSpace::FieldList<Dimension, Value> > {
  typedef FieldSpace::Field<Dimension, Value>* ElementType;
  static bool fixedSize() { return false; }
  static int numElements(const FieldSpace::FieldList<Dimension, Value>& x) { return x.size(); }
  static FieldSpace::Field<Dimension, Value> zero() { return FieldSpace::FieldList<Dimension, Value>(); }
};

}

#endif
