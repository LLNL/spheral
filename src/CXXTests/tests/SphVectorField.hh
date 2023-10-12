#ifndef SPHERAL_LVFIELD
#define SPHERAL_LVFIELD

#include "LvArray/bufferManipulation.hpp"
#include "Field/SphArray.hh"
#include "Field/FieldList.hh"


//template<typename Dimension, typename DATA_TYPE>
//class LvFieldView;

template<typename DataType>
class LvField;

template<typename DataType>
class SphVector:
  public chai::ManagedArray<DataType>{
  using ManagedArray = chai::ManagedArray<DataType>;

  friend LvField<DataType>;
public:

  RAJA_HOST SphVector() : ManagedArray(2) {}

  RAJA_HOST SphVector(size_t elems) : ManagedArray(elems), m_size(elems) {}

  RAJA_HOST void push_back(const DataType& value) {
    if (m_size >= capacity()) ManagedArray::reallocate(capacity() + (capacity() / 2));

    ManagedArray::data()[m_size] = value;
    m_size++;
  }

  RAJA_HOST void push_back(DataType&& value) {
    if (m_size >= capacity()) ManagedArray::reallocate(capacity() + (capacity() / 2));

    ManagedArray::data()[m_size] = std::move(value);
    m_size++;
  }

  template<typename... Args>
  RAJA_HOST
  DataType& emplace_back(Args&&... args) {
    if (m_size >= capacity()) ManagedArray::reallocate(capacity() + (capacity() / 2));

    new(&ManagedArray::data()[m_size]) DataType(std::forward<Args>(args)...);
    return ManagedArray::data()[m_size++];
  }
  
  RAJA_HOST_DEVICE size_t capacity() {return ManagedArray::m_elems;}
  RAJA_HOST_DEVICE size_t size() {return m_size;}

private:
  // *******************************************************
  // Required to Allow SphVector to be properly CHAICopyable
  RAJA_HOST_DEVICE SphVector<DataType>& operator=(std::nullptr_t) { ManagedArray::operator=(nullptr); return *this; }
  RAJA_HOST_DEVICE void shallowCopy(const SphVector& other) {
    m_size=other.m_size;
    ManagedArray::shallowCopy(other);
  }
  // *******************************************************

  size_t m_size = 0;
};

template<typename DATA_TYPE>
class LvField : public chai::CHAICopyable
{
public:

  using ARRAY_TYPE = Spheral::ManagedVector<DATA_TYPE>;

  RAJA_HOST LvField(size_t elems) : mDataArray{elems} {setName("");}

  RAJA_HOST LvField(size_t elems, RAJA::Platform platform) : mDataArray{elems} { 
    mDataArray.move(platform);
    setName("");
  }

  template<typename Dimension>
  RAJA_HOST LvField(std::string name, Spheral::NodeList<Dimension>& node_list) : mDataArray{node_list.numNodes()} { 
    setName(name);
  }

  RAJA_HOST LvField(std::string name, size_t elems, chai::ExecutionSpace const& platform = chai::CPU) : mDataArray{elems} { 
    strcpy(m_name, name.c_str());
    //mDataArray.move(platform);
    setName(name);
  }

  RAJA_HOST_DEVICE size_t size() const {return mDataArray.size();}

  template< typename U=LvField< DATA_TYPE > >
  RAJA_HOST
  void setName(std::string name) {
    strcpy(m_name, name.c_str());
    std::string const typeString = LvArray::system::demangleType< U >();
    mDataArray.setUserCallback(
      [typeString, name] (const chai::PointerRecord* record, chai::Action action, chai::ExecutionSpace exec) {
        if (action == chai::Action::ACTION_MOVE){
          std::string const size = LvArray::system::calculateSize(record->m_size);
          std::string const paddedSize = std::string( 9 - size.size(), ' ' ) + size;
          char const * const spaceStr = ( exec == chai::CPU ) ? "HOST  " : "DEVICE";
          std::cout << "Moved " << paddedSize << " to the " << spaceStr << ": " << typeString << " " << name <<std::endl;
        }
        if (action == chai::Action::ACTION_ALLOC){
          std::string const size = LvArray::system::calculateSize(record->m_size);
          std::string const paddedSize = std::string( 9 - size.size(), ' ' ) + size;
          std::cout << "Allocated " << paddedSize << " : " << typeString << " " << name <<std::endl;
        }
      }
    );
  }

  RAJA_HOST
  std::string name() const {return m_name;}

  RAJA_HOST
  //void move (RAJA::Platform platform) {mDataArray.move(platform);}
  void move(chai::ExecutionSpace const& space, bool touch = true) const {
    mDataArray.move(space,touch);
  }

  RAJA_HOST
  LvField make_pool_field(size_t num_pools, RAJA::Platform platform) const {
    return LvField(mDataArray.size() * num_pools, "POOL_" + name(), platform);
  }

  // Index operator.
  RAJA_HOST_DEVICE DATA_TYPE& operator[](const unsigned int index) {return mDataArray[index];}
  RAJA_HOST_DEVICE DATA_TYPE& operator[](const unsigned int index) const {return mDataArray[index];}

  RAJA_HOST_DEVICE inline constexpr LvField( LvField const & source) noexcept : mDataArray{source.mDataArray} {
  }
  ARRAY_TYPE& getArray() {return mDataArray;}

  // Required to Allow SphVector to be properly CHAICopyable
  RAJA_HOST_DEVICE LvField<DATA_TYPE>& operator=(std::nullptr_t) {mDataArray = nullptr; return *this;}
  RAJA_HOST_DEVICE void shallowCopy(const LvField& other) {
    strcpy(m_name, other.m_name);
    mDataArray.shallowCopy(other.mDataArray);
  }

private:
  ARRAY_TYPE mDataArray;
  char m_name[256] = "\0";

  friend LvField deepCopy(LvField const& rhs)
  {
    auto copy = LvField(
        rhs.m_name,
        rhs.size());
    copy.mDataArray = deepCopy(rhs.mDataArray);
    return copy;
  }
};


template<typename DIM, typename DATA_TYPE>
class LvFieldList{
public:


  //using ElementType = Spheral::Field<DIM,DATA_TYPE>;
  using ElementType = LvField<DATA_TYPE>;
  using ARRAY_TYPE = Spheral::ManagedVector<ElementType>;
  using ATOMIC_DATA_TYPE = typename DATA_TYPE::atomic_type;

  RAJA_HOST
  LvFieldList() {setName("\0");}

  RAJA_HOST
  LvFieldList(std::string name) {
    setName(name);
  }

  RAJA_HOST
  void appendField(const ElementType& field) { 
    mFieldArray.push_back(field);
  }

  RAJA_HOST
  void move(chai::ExecutionSpace const& space, bool touch = true) {
    mFieldArray.move(space, touch);
    for(int i = 0; i < mFieldArray.size(); i++) mFieldArray[i].move(space, touch);
  }

  RAJA_HOST_DEVICE
  DATA_TYPE& operator()(const unsigned int field_id, const unsigned int idx) const { return mFieldArray[field_id][idx]; }

  RAJA_HOST_DEVICE
  ATOMIC_DATA_TYPE& atomic(const unsigned field_id, const unsigned idx) const { return *reinterpret_cast<ATOMIC_DATA_TYPE*>(&mFieldArray[field_id][idx]); }

  RAJA_HOST_DEVICE
  ElementType& operator[](const unsigned int index) {return mFieldArray(index);}

  RAJA_HOST_DEVICE
  ElementType& operator[](const unsigned int index) const {return mFieldArray(index);}

  template< typename U=LvFieldList< DIM, DATA_TYPE > >
  RAJA_HOST
  void setName(std::string name) {
    strcpy(m_name, name.c_str());
    std::string const typeString = LvArray::system::demangleType< U >();
    mFieldArray.setUserCallback(
      [typeString, name] (const chai::PointerRecord* record, chai::Action action, chai::ExecutionSpace exec) {
        if (action == chai::Action::ACTION_MOVE){
          std::string const size = LvArray::system::calculateSize(record->m_size);
          std::string const paddedSize = std::string( 9 - size.size(), ' ' ) + size;
          char const * const spaceStr = ( exec == chai::CPU ) ? "HOST  " : "DEVICE";
          std::cout << "Moved " << paddedSize << " to the " << spaceStr << ": " << typeString << " " << name <<std::endl;
        }
      }
    );
  }

  RAJA_HOST_DEVICE
  inline LvFieldList( LvFieldList const & source) noexcept : mFieldArray{source.mFieldArray} {}

private:
  ARRAY_TYPE mFieldArray;
  char m_name[256] = "\0";
};


#endif

