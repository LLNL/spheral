#ifndef __Spheral_FieldView_hh__
#define __Spheral_FieldView_hh__

#include "config.hh"

namespace Spheral {

template <typename Dimension, typename DataType> class FieldView {
  DataType *m_data;
  size_t m_size;

public:
  SPHERAL_HOST_DEVICE
  FieldView(DataType *d, size_t s) : m_data(d), m_size(s) {}

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
  size_t size() const { return m_size; }

  SPHERAL_HOST_DEVICE
  DataType *data() const { return m_data; }
};

} // namespace Spheral

#endif
