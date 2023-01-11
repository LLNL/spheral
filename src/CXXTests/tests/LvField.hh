#ifndef SPHERAL_LVFIELD
#define SPHERAL_LVFIELD

#include "LvArray/Array.hpp"
#include "LvArray/ChaiBuffer.hpp"
#include "LvArray/bufferManipulation.hpp"

template<typename DATA_TYPE>
using LV_ARRAY_CHAI_1D = LvArray::Array<DATA_TYPE, 1, camp::idx_seq<0>, std::ptrdiff_t, LvArray::ChaiBuffer>;

template<typename DATA_TYPE>
using LV_ARRAY_VIEW_CHAI_1D = LvArray::ArrayView<DATA_TYPE, 1, 0, std::ptrdiff_t, LvArray::ChaiBuffer>;

template<typename DATA_TYPE>
class LvFieldView;

template<typename DATA_TYPE>
class LvFieldListView;


template<typename DATA_TYPE>
class LvField {
public:

  friend class LvFieldView<DATA_TYPE>;

  using ARRAY_TYPE = LV_ARRAY_CHAI_1D<DATA_TYPE>;
  using VIEW_TYPE = LvFieldView<DATA_TYPE>;

  LvField(size_t elems) { mDataArray = ARRAY_TYPE(elems); }

  LvField(size_t elems, RAJA::Platform platform) { 
    mDataArray = ARRAY_TYPE(elems);
    mDataArray.move(platform);
  }

  LvField(size_t elems, std::string name) { 
    mDataArray = ARRAY_TYPE(elems);
    mDataArray.setName(name);
  }

  LvField(size_t elems, std::string name, RAJA::Platform platform) { 
    mDataArray = ARRAY_TYPE(elems);
    mDataArray.setName(name);
    mDataArray.move(platform);
  }

  void setName(std::string name) {mDataArray.setName(name);}
  void move (RAJA::Platform platform) {mDataArray.move(platform);}
  LvFieldView<DATA_TYPE> toView() const {return LvFieldView<DATA_TYPE>(*this);}

  // Index operator.
  DATA_TYPE& operator[](const unsigned int index) {return mDataArray[index];}
  DATA_TYPE& operator[](const unsigned int index) const {return mDataArray[index];}

  ARRAY_TYPE& getArray() {return mDataArray;}

private:
  ARRAY_TYPE mDataArray;

};

template<typename DATA_TYPE>
class LvFieldView {
public:
  using ARRAY_VIEW_TYPE = LV_ARRAY_VIEW_CHAI_1D<DATA_TYPE>;
  using FIELD_TYPE = LvField<DATA_TYPE>;

  friend class LvFieldListView<DATA_TYPE>;

  LvFieldView(const FIELD_TYPE& field) : mDataView{field.mDataArray} {}

  void move(LvArray::MemorySpace const& space, bool touch = true) const {
    mDataView.move(space,touch);
  }

  ARRAY_VIEW_TYPE& getView() {return mDataView;}

  RAJA_HOST_DEVICE
  DATA_TYPE& operator[](const unsigned int index) {return mDataView(index);}

  RAJA_HOST_DEVICE
  DATA_TYPE& operator[](const unsigned int index) const {return mDataView(index);}

  RAJA_HOST_DEVICE
  inline constexpr LvFieldView( LvFieldView const & source) noexcept : mDataView{source.mDataView} {}

private:
  ARRAY_VIEW_TYPE mDataView;
};


template<typename DATA_TYPE>
class LvFieldList{
public:

  friend class LvFieldListView<DATA_TYPE>;

  using FIELD_TYPE = LvField<DATA_TYPE>;
  using FIELD_VIEW_TYPE = LvFieldView<DATA_TYPE>;
  using ARRAY_TYPE = LV_ARRAY_CHAI_1D<FIELD_VIEW_TYPE>;

  LvFieldList() {}

  LvFieldList(std::string name) {
    mFieldArray.setName(name);
  }

  void appendField(const FIELD_TYPE& field) { 
    mFieldArray.emplace_back(field.toView());
  }

  void appendField(const FIELD_VIEW_TYPE& field) { 
    mFieldArray.emplace_back(field);
  }

private:
  ARRAY_TYPE mFieldArray;
};


template<typename DATA_TYPE>
class LvFieldListView{
public:

  using FIELD_VIEW_TYPE = LvFieldView<DATA_TYPE>;
  using ARRAY_VIEW_TYPE = LV_ARRAY_VIEW_CHAI_1D<FIELD_VIEW_TYPE>;

  RAJA_HOST_DEVICE
  LvFieldListView(const LvFieldList<DATA_TYPE>& field) : mFieldView{field.mFieldArray.toView()} {}

  void move(LvArray::MemorySpace const& space, bool touch = true) const {
    mFieldView.move(space,touch);
  }

  RAJA_HOST_DEVICE
  DATA_TYPE& operator()(const unsigned int field_id, const unsigned int idx) const { return mFieldView(field_id)[idx]; }

  RAJA_HOST_DEVICE
  FIELD_VIEW_TYPE& operator[](const unsigned int index) {return mFieldView(index);}

  RAJA_HOST_DEVICE
  FIELD_VIEW_TYPE& operator[](const unsigned int index) const {return mFieldView(index);}

  RAJA_HOST_DEVICE
  inline LvFieldListView( LvFieldListView const & source) noexcept : mFieldView{source.mFieldView} {}

private:
  ARRAY_VIEW_TYPE mFieldView;
};


#endif

