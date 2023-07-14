#ifndef __Spheral_Field_View_hh__
#define __Spheral_Field_View_hh__

#include "Field/SphArray.hh"
#include "Field.hh"

template<typename Dimension, typename DataType>
class FieldListView;

namespace Spheral {

template<typename Dimension, typename DataType>
class FieldView {
public:
  using ArrayViewType = SphArrayView<DataType>;
  using FieldType = Field<Dimension, DataType>;

  friend class FieldListView<Dimension, DataType>;

  FieldView();
  FieldView(const FieldType& field);
  FieldView(const FieldType& field, const FieldType& pool);

  void move(LvArray::MemorySpace const& space, bool touch = true) const;

  ArrayViewType& getView();
  FieldType& get() const;
  FieldType& operator*() const;
  FieldType* operator->() const;

  RAJA_HOST_DEVICE
  DataType& operator[](const unsigned int index);

  RAJA_HOST_DEVICE
  DataType& operator[](const unsigned int index) const;

  RAJA_HOST_DEVICE
  DataType& pool(const unsigned index) const;

  RAJA_HOST_DEVICE
  inline constexpr FieldView( FieldView const & source) noexcept;

  FieldView& operator=(const FieldView& rhs);
  bool operator==(const FieldType& rhs) const;
  bool operator==(const FieldView& rhs) const;

private:
  ArrayViewType mDataView;
  ArrayViewType mDataPoolView;
  FieldType* mFieldPtr = nullptr;
};

} // namespace Spheral

#include "FieldViewInline.hh"

#else //  __Spheral_Field_View_hh__

// Forward declare the Field class.
namespace Spheral {
  template<typename Dimension, typename DataType> class FieldView;
} // namespace Spheral


#endif //  __Spheral_Field_View_hh__


