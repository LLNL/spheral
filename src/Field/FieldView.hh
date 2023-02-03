#ifndef __Spheral_Field_View_hh__
#define __Spheral_Field_View_hh__

#include "Field/SphArray.hh"
#include "Field.hh"

template<typename Dimension, typename DataType>
class LvFieldListView;

namespace Spheral {

template<typename Dimension, typename DataType>
class FieldView {
public:
  using ARRAY_VIEW_TYPE = SphArrayView<DataType>;
  using FIELD_TYPE = Field<Dimension, DataType>;

  friend class ::LvFieldListView<Dimension, DataType>;

  FieldView();
  FieldView(const FIELD_TYPE& field);
  FieldView(const FIELD_TYPE& field, const FIELD_TYPE& pool);

  void move(LvArray::MemorySpace const& space, bool touch = true) const;

  ARRAY_VIEW_TYPE& getView();
  FIELD_TYPE& get() const;
  FIELD_TYPE& operator*() const;
  FIELD_TYPE* operator->() const;

  RAJA_HOST_DEVICE
  DataType& operator[](const unsigned int index);

  RAJA_HOST_DEVICE
  DataType& operator[](const unsigned int index) const;

  RAJA_HOST_DEVICE
  DataType& pool(const unsigned index) const;

  RAJA_HOST_DEVICE
  inline constexpr FieldView( FieldView const & source) noexcept;

  FieldView& operator=(const FieldView& rhs);
  bool operator==(const FIELD_TYPE& rhs) const;
  bool operator==(const FieldView& rhs) const;

private:
  ARRAY_VIEW_TYPE mDataView;
  ARRAY_VIEW_TYPE mDataPoolView;
  FIELD_TYPE* mFieldPtr = nullptr;
};

} // namespace Spheral

#include "FieldViewInline.hh"

#else //  __Spheral_Field_View_hh__

// Forward declare the Field class.
namespace Spheral {
  template<typename Dimension, typename DataType> class FieldView;
} // namespace Spheral


#endif //  __Spheral_Field_View_hh__


