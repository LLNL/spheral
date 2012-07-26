//---------------------------------Spheral++----------------------------------//
// FieldListSet -- A container class for FieldLists.  Provides a convenient
//   container for a set of associated FieldLists, allowing us to pass them
//   in bundles through routines such as FieldListFunctions.
//
// Created by JMO, Wed Jun  5 21:47:20 PDT 2002
//----------------------------------------------------------------------------//
#ifndef __Spheral_FieldListSet_hh__
#define __Spheral_FieldListSet_hh__

#ifndef __GCCXML__
#include <vector>
#include "FieldList.hh"
#else
#include "fakestl.hh"
#endif

namespace Spheral {
namespace FieldSpace {

template<typename Dimension>
struct FieldListSet {
  std::vector< FieldList<Dimension, typename Dimension::Scalar> > ScalarFieldLists;
  std::vector< FieldList<Dimension, typename Dimension::Vector> > VectorFieldLists;
  std::vector< FieldList<Dimension, typename Dimension::Tensor> > TensorFieldLists;
  std::vector< FieldList<Dimension, typename Dimension::SymTensor> > SymTensorFieldLists;

  // These methods are here as helpers for pybindgen.
  std::vector< FieldList<Dimension, typename Dimension::Scalar> >* ScalarFieldListPtrs();
  std::vector< FieldList<Dimension, typename Dimension::Vector> >* VectorFieldListPtrs();
  std::vector< FieldList<Dimension, typename Dimension::Tensor> >* TensorFieldListPtrs();
  std::vector< FieldList<Dimension, typename Dimension::SymTensor> >* SymTensorFieldListPtrs();
};

}
}

#ifndef __GCCXML__
#include "FieldListSetInline.hh"
#endif

#else

namespace Spheral {
  namespace FieldSpace {
    // Forward declaration.
    template<typename Dimension> struct FieldListSet;
  }
}

#endif
