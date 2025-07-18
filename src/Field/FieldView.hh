#ifndef __Spheral_FieldView_hh__
#define __Spheral_FieldView_hh__

#include "chai/ManagedArray.hpp"
#include "config.hh"

namespace Spheral {

template <typename Dimension, typename DataType> class FieldView {


  using ContainerType = typename chai::ManagedArray<DataType>;

  ContainerType mData;

public:
  SPHERAL_HOST_DEVICE FieldView() {}

  SPHERAL_HOST
  FieldView(ContainerType const &d) : mData(d) {}

  SPHERAL_HOST_DEVICE
  FieldView(FieldView const &) = default;
  SPHERAL_HOST_DEVICE
  FieldView(FieldView &&) = default;
  SPHERAL_HOST_DEVICE
  FieldView& operator=(const FieldView&) = default;

  // Element access.
  SPHERAL_HOST_DEVICE
  DataType &operator()(int index) { return FieldView::operator[](index); }

  SPHERAL_HOST_DEVICE
  DataType &at(int index) { return FieldView::operator[](index); }

  // Index operator.
  SPHERAL_HOST_DEVICE
  DataType &operator[](const size_t index) { return mData[index]; }
  SPHERAL_HOST_DEVICE
  DataType &operator[](const size_t index) const { return mData[index]; }

  SPHERAL_HOST_DEVICE
  size_t size() const { return mData.size(); }

  SPHERAL_HOST_DEVICE
  DataType *data() const { return mData.getActivePointer(); }
};

} // namespace Spheral

#endif
