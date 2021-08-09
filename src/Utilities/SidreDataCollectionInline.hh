//SidreDataCollectionInLine.hh
#include "axom/sidre.hpp"
#include "Field/Field.hh"

namespace Spheral
{

template <typename Dimension, typename DataType,
         typename std::enable_if<std::is_arithmetic<DataType>::value,
                                 DataType>::type* = nullptr>
inline
axom::sidre::Group *SidreDataCollection::sidreStoreField(const std::string &view_name, 
                                                         const Spheral::Field<Dimension, DataType> &field)
{
   axom::sidre::DataTypeId dtype = field.getAxomTypeID();
   axom::sidre::Group* wholeField = m_datastore_ptr->getRoot()->createGroup(view_name + "_Group");
   axom::IndexType num_elements = 1;
   int view_count = 0;

   for (u_int i = 0; i < field.size(); ++i)
   {
      const DataType *data = &(field[i]);
      wholeField->createView(view_name + std::to_string(view_count), dtype, num_elements, (void*)data);
      view_count++;
   }
   return wholeField;
}

//------------------------------------------------------------------------------
template<typename Dimension>
inline
axom::sidre::Group *SidreDataCollection::sidreStoreField(const std::string &view_name, 
                                                         const Spheral::Field<Dimension, std::string> &field)
{
   axom::sidre::DataTypeId dtype = field.getAxomTypeID();
   axom::sidre::Group* wholeField = m_datastore_ptr->getRoot()->createGroup(view_name + "_Group");
   int view_count = 0;

   for (u_int i = 0; i < field.size(); ++i)
   {
      axom::IndexType num_elements = field[i].size();
      const char *data = &(*field[i].begin());
      wholeField->createView(view_name + std::to_string(view_count), dtype, num_elements, (void*)data);
      view_count++;
   }
   return wholeField;
}

//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
axom::sidre::Group *SidreDataCollection::sidreStoreField(const std::string &view_name, 
                                                         const Spheral::Field<Dimension, std::vector<DataType>> &field)
{
   axom::sidre::DataTypeId dtype = field.getAxomTypeID();
   axom::sidre::Group* wholeField = m_datastore_ptr->getRoot()->createGroup(view_name + "_Group");
   int view_count = 0;

   for (u_int i = 0; i < field.size(); ++i)
   {
      axom::IndexType num_elements = field[i].size();
      const DataType *data = &(*field[i].begin());
      wholeField->createView(view_name + std::to_string(view_count), dtype, num_elements, (void*)data);
      view_count++;
   }
   return wholeField;
}

//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
axom::sidre::Group *SidreDataCollection::sidreStoreField(const std::string &view_name,
                                                         const Spheral::Field<Dimension, std::tuple<DataType, DataType, DataType>> &field)
{
   axom::sidre::DataTypeId dtype = field.getAxomTypeID();
   axom::sidre::Group* wholeField = m_datastore_ptr->getRoot()->createGroup(view_name + "_Group");
   axom::IndexType num_elements = 3;
   int view_count = 0;

   for (u_int i = 0; i < field.size(); ++i)
   {
      DataType data [] = {std::get<0>(field[i]), std::get<1>(field[i]), std::get<2>(field[i])};
      axom::sidre::Buffer* buff = m_datastore_ptr->createBuffer()->allocate(dtype, num_elements)
                                                 ->copyBytesIntoBuffer(data, sizeof(DataType) * num_elements);
      wholeField->createView(view_name + std::to_string(view_count), dtype, num_elements, buff);
      view_count++;
   }
   return wholeField;
}

//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
axom::sidre::Group *SidreDataCollection::sidreStoreField(const std::string &view_name,
                                                         const Spheral::Field<Dimension, std::tuple<DataType, DataType, DataType, DataType>> &field)
{
   axom::sidre::DataTypeId dtype = field.getAxomTypeID();
   axom::sidre::Group* wholeField = m_datastore_ptr->getRoot()->createGroup(view_name + "_Group");
   axom::IndexType num_elements = 4;
   int view_count = 0;

   for (u_int i = 0; i < field.size(); ++i)
   {
      DataType data [] = {std::get<0>(field[i]), std::get<1>(field[i]), std::get<2>(field[i]), std::get<3>(field[i])};
      axom::sidre::Buffer* buff = m_datastore_ptr->createBuffer()->allocate(dtype, num_elements)
                                                 ->copyBytesIntoBuffer(data, sizeof(DataType) * num_elements);
      wholeField->createView(view_name + std::to_string(view_count), dtype, num_elements, buff);
      view_count++;
   }
   return wholeField;
}

//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
axom::sidre::Group *SidreDataCollection::sidreStoreField(const std::string &view_name,
                                                         const Spheral::Field<Dimension, std::tuple<DataType, DataType, DataType, DataType, DataType>> &field)
{
   axom::sidre::DataTypeId dtype = field.getAxomTypeID();
   axom::sidre::Group* wholeField = m_datastore_ptr->getRoot()->createGroup(view_name + "_Group");
   axom::IndexType num_elements = 5;
   int view_count = 0;

   for (u_int i = 0; i < field.size(); ++i)
   {
      DataType data [] = {std::get<0>(field[i]), std::get<1>(field[i]), std::get<2>(field[i]), std::get<3>(field[i]), std::get<4>(field[i])};
      axom::sidre::Buffer* buff = m_datastore_ptr->createBuffer()->allocate(dtype, num_elements)
                                                 ->copyBytesIntoBuffer(data, sizeof(DataType) * num_elements);
      wholeField->createView(view_name + std::to_string(view_count), dtype, num_elements, buff);
      view_count++;
   }
   return wholeField;
}

//------------------------------------------------------------------------------
template <typename Dimension, typename DataType,
         typename std::enable_if<!is_rank_n_tensor<DataType>::value && !std::is_arithmetic<DataType>::value,
                                 DataType>::type* = nullptr>
inline
axom::sidre::Group *SidreDataCollection::sidreStoreField(const std::string &view_name, 
                                                         const Spheral::Field<Dimension, DataType> &field)
{
   axom::sidre::DataTypeId dtype = field.getAxomTypeID();
   axom::sidre::Group* wholeField = m_datastore_ptr->getRoot()->createGroup(view_name + "_Group");
   axom::IndexType num_elements = DataTypeTraits<DataType>::numElements(field[0]);
   int view_count = 0;

   for (u_int i = 0; i < field.size(); ++i)
   {
      auto *data = &(*field[i].begin());
      wholeField->createView(view_name + std::to_string(view_count), dtype, num_elements, (void*)data);
      view_count++;
   }
   return wholeField;
}

template <typename Dimension, typename DataType,
         typename std::enable_if<is_rank_n_tensor<DataType>::value && !std::is_arithmetic<DataType>::value,
                                 DataType>::type* = nullptr>
inline
axom::sidre::Group *SidreDataCollection::sidreStoreField(const std::string &view_name, 
                                                         const Spheral::Field<Dimension, DataType> &field)
{
   axom::sidre::DataTypeId dtype = field.getAxomTypeID();
   axom::sidre::Group* wholeField = m_datastore_ptr->getRoot()->createGroup(view_name + "_Group");
   axom::IndexType num_elements = DataTypeTraits<DataType>::numElements(field[0]);
   int view_count = 0;

   double data [num_elements];
   for (u_int i = 0; i < field.size(); ++i)
   {
      int index = 0;
      for (auto it = field[i].begin(); it != field[i].end(); ++it)
      {
         data[index] = *it;
         index++;
      }
      axom::sidre::Buffer* buff = m_datastore_ptr->createBuffer()->allocate(dtype, num_elements)
                                                 ->copyBytesIntoBuffer(data, sizeof(double) * num_elements);
      wholeField->createView(view_name + std::to_string(view_count), dtype, num_elements, buff);
      view_count++;
   }
   return wholeField;
}


template<typename Dimension>
inline
axom::sidre::Group *SidreDataCollection::sidreStoreField(const std::string &view_name, 
                                                         const Spheral::Field<Dimension, Dim<2>::Vector> &field)
{
   axom::sidre::DataTypeId dtype = field.getAxomTypeID();
   axom::sidre::Group* wholeField = m_datastore_ptr->getRoot()->createGroup(view_name + "_Group");
   axom::IndexType num_elements = DataTypeTraits<Dim<2>::Vector>::numElements(field[0]);
   int view_count = 0;

   double data [num_elements];
   for (u_int i = 0; i < field.size(); ++i)
   {
      int index = 0;
      for (auto it = field[i].begin(); it != field[i].end(); ++it)
      {
         data[index] = *it;
         index++;
      }
      axom::sidre::Buffer* buff = m_datastore_ptr->createBuffer()->allocate(dtype, num_elements)
                                                 ->copyBytesIntoBuffer(data, sizeof(double) * num_elements);
      wholeField->createView(view_name + std::to_string(view_count), dtype, num_elements, buff);
      view_count++;
   }
   return wholeField;
}

} // namespace Spheral
