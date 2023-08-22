#ifndef SPHERAL_LVFIELD
#define SPHERAL_LVFIELD

#include "LvArray/bufferManipulation.hpp"
#include "Field/SphArray.hh"
#include "Field/FieldView.hh"
#include "Field/FieldList.hh"
#include "Field/FieldListView.hh"



//template<typename Dimension, typename DATA_TYPE>
//class LvFieldView;

//template<typename DATA_TYPE>
//class LvField {
//public:
//
//  friend class LvFieldView<DATA_TYPE>;
//
//  using ARRAY_TYPE = Spheral::SphArray<DATA_TYPE>;
//  using VIEW_TYPE = LvFieldView<DATA_TYPE>;
//
//  LvField(size_t elems) { mDataArray = ARRAY_TYPE(elems); }
//
//  LvField(size_t elems, RAJA::Platform platform) { 
//    mDataArray = ARRAY_TYPE(elems);
//    mDataArray.move(platform);
//  }
//
//  LvField(size_t elems, std::string name) { 
//    mDataArray = ARRAY_TYPE(elems);
//    mName = name; 
//    mDataArray.setName(mName);
//  }
//
//  LvField(size_t elems, std::string name, RAJA::Platform platform) { 
//    mDataArray = ARRAY_TYPE(elems);
//    mName = name; 
//    mDataArray.setName(mName);
//    mDataArray.move(platform);
//  }
//
//  void setName(std::string name) {mDataArray.setName(name);}
//  void move (RAJA::Platform platform) {mDataArray.move(platform);}
//
//  LvFieldView<DATA_TYPE> toView() const {return LvFieldView<DATA_TYPE>(*this);}
//  LvFieldView<DATA_TYPE> toViewWithPool(const LvField& pool) const {return LvFieldView<DATA_TYPE>(*this, pool);}
//
//  LvField make_pool_field(size_t num_pools, RAJA::Platform platform) {
//    return LvField(mDataArray.size() * num_pools, "POOL_" + mName, platform);
//  }
//
//  // Index operator.
//  DATA_TYPE& operator[](const unsigned int index) {return mDataArray[index];}
//  DATA_TYPE& operator[](const unsigned int index) const {return mDataArray[index];}
//
//  ARRAY_TYPE& getArray() {return mDataArray;}
//
//private:
//  ARRAY_TYPE mDataArray;
//  std::string mName = "";
//
//};

//template<typename Dimension, typename DATA_TYPE>
//class LvFieldView {
//public:
//  using ARRAY_VIEW_TYPE = Spheral::SphArrayView<DATA_TYPE>;
//  using FIELD_TYPE = Spheral::Field<Dimension, DATA_TYPE>;
//
//  friend class LvFieldListView<Dimension, DATA_TYPE>;
//
//  LvFieldView(const FIELD_TYPE& field) : mDataView{field.mDataArray}, mDataPoolView{field.mDataArray} {}
//  LvFieldView(const FIELD_TYPE& field, const FIELD_TYPE& pool) : mDataView{field.mDataArray}, mDataPoolView{pool.mDataArray} {}
//
//  void move(LvArray::MemorySpace const& space, bool touch = true) const {
//    mDataView.move(space,touch);
//  }
//
//  ARRAY_VIEW_TYPE& getView() {return mDataView;}
//
//  RAJA_HOST_DEVICE
//  DATA_TYPE& operator[](const unsigned int index) {return mDataView(index); }
//
//  RAJA_HOST_DEVICE
//  DATA_TYPE& operator[](const unsigned int index) const {return mDataView(index); }
//
//  RAJA_HOST_DEVICE
//  DATA_TYPE& pool(const unsigned index) const {return mDataPoolView(index); }
//
//  RAJA_HOST_DEVICE
//  inline constexpr LvFieldView( LvFieldView const & source) noexcept : mDataView{source.mDataView}, mDataPoolView{source.mDataPoolView} {}
//
//private:
//  ARRAY_VIEW_TYPE mDataView;
//  ARRAY_VIEW_TYPE mDataPoolView;
//};


//template<typename Dimension, typename DATA_TYPE>
//class LvFieldList{
//public:
//
//  friend class LvFieldListView<Dimension, DATA_TYPE>;
//
//  using FIELD_TYPE = Spheral::Field<Dimension, DATA_TYPE>;
//  //using FIELD_TYPE = LvField<DATA_TYPE>;
//  using FIELD_VIEW_TYPE = Spheral::FieldView<Dimension, DATA_TYPE>;
//  using ARRAY_TYPE = Spheral::SphArray<FIELD_VIEW_TYPE>;
//
//  LvFieldList() {}
//
//  LvFieldList(std::string name) {
//    mFieldArray.setName(name);
//  }
//
//  void appendField(const FIELD_TYPE& field) { 
//    mFieldArray.emplace_back(field.toView());
//  }
//
//  void appendField(const FIELD_TYPE& field, const FIELD_TYPE& pool) { 
//    mFieldArray.emplace_back(field.toViewWithPool(pool));
//  }
//
//  void appendField(const FIELD_VIEW_TYPE& field) { 
//    mFieldArray.emplace_back(field);
//  }
//
//private:
//  ARRAY_TYPE mFieldArray;
//};


//template<typename Dimension, typename DATA_TYPE>
//class LvFieldListView{
//public:
//
//  using FIELD_VIEW_TYPE = Spheral::FieldView<Dimension, DATA_TYPE>;
//  using ARRAY_VIEW_TYPE = Spheral::SphArrayView<FIELD_VIEW_TYPE>;
//
//#if USE_DEVICE
//  using ATOMIC_DATA_TYPE = typename DATA_TYPE::atomic_type;
//#else
//  using ATOMIC_DATA_TYPE = DATA_TYPE;
//#endif
//
//  RAJA_HOST_DEVICE
//  LvFieldListView(const Spheral::FieldList<Dimension, DATA_TYPE>& field) : mFieldViewsView{field.mFieldPtrs.toView()} {}
//
//  void move(LvArray::MemorySpace const& space, bool touch = true) const {
//    mFieldViewsView.move(space,touch);
//  }
//
//  RAJA_HOST_DEVICE
//  DATA_TYPE& operator()(const unsigned int field_id, const unsigned int idx) const { return mFieldViewsView(field_id)[idx]; }
//
//  RAJA_HOST_DEVICE
//  FIELD_VIEW_TYPE& operator[](const unsigned int index) {return mFieldViewsView(index);}
//
//  RAJA_HOST_DEVICE
//  FIELD_VIEW_TYPE& operator[](const unsigned int index) const {return mFieldViewsView(index);}
//
//  RAJA_HOST_DEVICE
//  DATA_TYPE& pool(const unsigned field_id, const unsigned idx) const {return mFieldViewsView(field_id).pool(idx); }
//
//  RAJA_HOST_DEVICE
//  ATOMIC_DATA_TYPE& pool_atomic(const unsigned field_id, const unsigned idx) const { return *reinterpret_cast<ATOMIC_DATA_TYPE*>(&mFieldViewsView(field_id).pool(idx)); }
//
//  RAJA_HOST_DEVICE
//  inline LvFieldListView( LvFieldListView const & source) noexcept : mFieldViewsView{source.mFieldViewsView} {}
//
//private:
//  ARRAY_VIEW_TYPE mFieldViewsView;
//};


#endif

