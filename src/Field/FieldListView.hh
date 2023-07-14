#ifndef __Spheral_FieldList_View_hh__
#define __Spheral_FieldList_View_hh__

#include "Field/SphArray.hh"
#include "FieldList.hh"

//template<typename Dimension, typename DataType>
//class FieldListView;

namespace Spheral {

template<typename Dimension, typename DataType>
class FieldListView{
public:

  using FieldViewType = FieldView<Dimension, DataType>;
  using ArrayViewType = SphArrayView<FieldViewType>;

#if USE_DEVICE
  using AtomicDataType = typename DataType::atomic_type;
#else
  using AtomicDataType = DataType;
#endif

  RAJA_HOST_DEVICE
  FieldListView(const FieldList<Dimension, DataType>& field);

  void move(LvArray::MemorySpace const& space, bool touch = true) const;

  RAJA_HOST_DEVICE
  DataType& operator()(const unsigned int field_id, const unsigned int idx) const;

  RAJA_HOST_DEVICE
  FieldViewType& operator[](const unsigned int index);

  RAJA_HOST_DEVICE
  FieldViewType& operator[](const unsigned int index) const;

  RAJA_HOST_DEVICE
  DataType& pool(const unsigned field_id, const unsigned idx) const;

  RAJA_HOST_DEVICE
  AtomicDataType& atomic(const unsigned field_id, const unsigned idx) const;

  RAJA_HOST_DEVICE
  AtomicDataType& pool_atomic(const unsigned field_id, const unsigned idx) const;

  RAJA_HOST_DEVICE
  inline FieldListView( FieldListView const & source) noexcept;

private:
  ArrayViewType mFieldViewsView;
};


} // namespace Spheral

#include "FieldListViewInline.hh"

#else //  __Spheral_Field_View_hh__

// Forward declare the Field class.
namespace Spheral {
  template<typename Dimension, typename DataType> class FieldListView;
} // namespace Spheral
  //
#endif
