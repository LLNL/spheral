#ifndef __PBGWrap_FieldOperationsType__
#define __PBGWrap_FieldOperationsType__

#include "FieldOperations/FieldListFunctions.hh"
#include "FieldOperations/FieldListFunctionsMash.hh"
#include "FieldOperations/FieldListSecondDerivatives.hh"
#include "FieldOperations/PairWiseFieldListFunctions.hh"
#include "FieldOperations/sampleMultipleFields2Lattice.hh"
#include "FieldOperations/binFieldList2Lattice.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Disambiguate the binFields methods.
//------------------------------------------------------------------------------
template<typename Dimension, typename Value>
inline
std::vector<Value>
binFieldList2LatticeWithSmoothing(const FieldList<Dimension, Value>& fieldList,
                                  const TableKernel<Dimension>& W,
                                  const typename Dimension::Vector& xmin,
                                  const typename Dimension::Vector& xmax,
                                  const std::vector<unsigned>& nsample) {
  return binFieldList2Lattice(fieldList, W, xmin, xmax, nsample);
}

}

#endif
