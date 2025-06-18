#ifndef __Spheral_FieldView_hh__
#define __Spheral_FieldView_hh__

#include "chai/ManagedArray.hpp"
#include "config.hh"

namespace Spheral {

template <typename Dimension, typename DataType> class FieldView {

  // Ensure the datatype is trivially copyable.
  static_assert(std::is_trivially_copyable<DataType>::value,
                "Error: The template type T must be trivially copyable.");

  using ContainerType = typename chai::ManagedArray<DataType>;

  ContainerType m_data = chai::ArrayManager::s_null_record;

public:
  SPHERAL_HOST_DEVICE FieldView() {}

  SPHERAL_HOST
  FieldView(ContainerType const &d) : m_data(d) {}

  SPHERAL_HOST_DEVICE
  FieldView(FieldView const &) = default;
  SPHERAL_HOST_DEVICE
  FieldView(FieldView &&) = default;

  // Element access.
  SPHERAL_HOST_DEVICE
  DataType &operator()(int index) { return FieldView::operator[](index); }

  SPHERAL_HOST_DEVICE
  DataType &at(int index) { return FieldView::operator[](index); }

  // Index operator.
  SPHERAL_HOST_DEVICE
  DataType &operator[](const size_t index) { return m_data[index]; }
  SPHERAL_HOST_DEVICE
  DataType &operator[](const size_t index) const { return m_data[index]; }

  SPHERAL_HOST_DEVICE
  size_t size() const { return m_data.size(); }

  SPHERAL_HOST_DEVICE
  DataType *data() const { return m_data.getActivePointer(); }
};

} // namespace Spheral

#endif
