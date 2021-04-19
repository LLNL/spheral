#include "axom/sidre.hpp"
#include "Field/Field.hh"

namespace Spheral
{

template<typename Dimension, typename DataType>
inline
axom::sidre::View *SidreDataCollection::alloc_view(const std::string &view_name, 
                                                   const Spheral::Field<Dimension, DataType> &field)
{
   axom::sidre::DataTypeId dtype = field.getAxomType();
   axom::IndexType num_elements = field.numElements();
   const DataType *data = &(*field.begin());
   axom::sidre::View *v = m_datastore_ptr->getRoot()->createView(view_name, dtype,
                                                        num_elements, (void*)data);
   return v;
}

//------------------------------------------------------------------------------
template<typename Dimension>
inline
axom::sidre::View *SidreDataCollection::alloc_view(const std::string &view_name, 
                                                   const Spheral::Field<Dimension, std::string> &field)
{
   axom::sidre::DataTypeId dtype = field.getAxomType();
   int view_count = 0;
   for (u_int i = 0; i < field.size(); i++)
   {
      axom::IndexType num_elements = field[i].size();
      const char *data = &(*field[i].begin());
      m_datastore_ptr->getRoot()->createView(view_name + std::to_string(view_count), dtype, 
                                                         num_elements, (void*)data);
      view_count++;
   }
   return m_datastore_ptr->getRoot()->getView(view_name + "0");
}

//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
axom::sidre::View *SidreDataCollection::alloc_view(const std::string &view_name, 
                                                   const Spheral::Field<Dimension, std::vector<DataType>> &field)
{
   axom::sidre::DataTypeId dtype = field.getAxomType();
   int view_count = 0;
   for (u_int i = 0; i < field.size(); i++)
   {
      axom::IndexType num_elements = field[i].size();
      const DataType *data = &(*field[i].begin());
      m_datastore_ptr->getRoot()->createView(view_name + std::to_string(view_count), dtype, 
                                                         num_elements, (void*)data);
      view_count++;
   }
   return m_datastore_ptr->getRoot()->getView(view_name + "0");
}

//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
axom::sidre::View *SidreDataCollection::alloc_view(const std::string &view_name,
                              const Spheral::Field<Dimension, std::tuple<DataType, DataType, DataType>> &field)
{
   axom::sidre::DataTypeId dtype = field.getAxomType();
   axom::IndexType num_elements = 3; //field.numElements(field[0]); //std::tuple_size
   // std::cout << field.numElements() << std::endl;
   int view_count = 0;

   for (u_int i = 0; i < field.size(); i++)
   {
      DataType data [] = {std::get<0>(field[i]), std::get<1>(field[i]), std::get<2>(field[i])};
      axom::sidre::Buffer* buff = m_datastore_ptr->createBuffer()->allocate(dtype, num_elements)
                                                 ->copyBytesIntoBuffer(data, sizeof(DataType) * num_elements);
      m_datastore_ptr->getRoot()->createView(view_name + std::to_string(view_count), dtype, 
                                                         num_elements, buff);
      view_count++;
   }
   return m_datastore_ptr->getRoot()->getView(view_name + "0");
}

//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
axom::sidre::View *SidreDataCollection::alloc_view(const std::string &view_name,
                              const Spheral::Field<Dimension, std::tuple<DataType, DataType, DataType, DataType>> &field)
{
   axom::sidre::DataTypeId dtype = field.getAxomType();
   axom::IndexType num_elements = 4;
   int view_count = 0;
   for (u_int i = 0; i < field.size(); i++)
   {
      DataType data [] = {std::get<0>(field[i]), std::get<1>(field[i]), std::get<2>(field[i]), std::get<3>(field[i])};
      axom::sidre::Buffer* buff = m_datastore_ptr->createBuffer()->allocate(dtype, num_elements)
                                                 ->copyBytesIntoBuffer(data, sizeof(DataType) * num_elements);
      m_datastore_ptr->getRoot()->createView(view_name + std::to_string(view_count), dtype, 
                                                         num_elements, buff);
      view_count++;
   }
   return m_datastore_ptr->getRoot()->getView(view_name + "0");
}

//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
axom::sidre::View *SidreDataCollection::alloc_view(const std::string &view_name,
                              const Spheral::Field<Dimension, std::tuple<DataType, DataType, DataType, DataType, DataType>> &field)
{
   axom::sidre::DataTypeId dtype = field.getAxomType();
   axom::IndexType num_elements = 5;
   int view_count = 0;
   for (u_int i = 0; i < field.size(); i++)
   {
      DataType data [] = {std::get<0>(field[i]), std::get<1>(field[i]), std::get<2>(field[i]), std::get<3>(field[i]), std::get<4>(field[i])};
      axom::sidre::Buffer* buff = m_datastore_ptr->createBuffer()->allocate(dtype, num_elements)
                                                 ->copyBytesIntoBuffer(data, sizeof(DataType) * num_elements);
      m_datastore_ptr->getRoot()->createView(view_name + std::to_string(view_count), dtype, 
                                                         num_elements, buff);
      view_count++;
   }
   return m_datastore_ptr->getRoot()->getView(view_name + "0");
}

// //------------------------------------------------------------------------------
// template<typename Dimension>
// inline
// axom::sidre::View *SidreDataCollection::alloc_view(const std::string &view_name, 
//                                                    const Spheral::Field<Dimension, std::pair<Value1, Value2>> &field)
// {
//    axom::sidre::DataTypeId dtype = field.getAxomType();
//    int view_count = 0;
//    for (u_int i = 0; i < field.size(); i++)
//    {
//       axom::IndexType num_elements = field[i].size();
//       const char *data = &(*field[i].begin());
//       m_datastore_ptr->getRoot()->createView(view_name + std::to_string(view_count), dtype, 
//                                                          num_elements, (void*)data);
//       view_count++;
//    }
//    return m_datastore_ptr->getRoot()->getView(view_name + "0");
// }

// //------------------------------------------------------------------------------
// template<typename Dimension>
// inline
// axom::sidre::View *SidreDataCollection::alloc_view(const std::string &view_name, 
//                                                    const Spheral::Field<Dimension, std::unordered_map<Value1, Value2, Hash>> &field)
// {
//    axom::sidre::DataTypeId dtype = field.getAxomType();
//    int view_count = 0;
//    for (u_int i = 0; i < field.size(); i++)
//    {
//       axom::IndexType num_elements = field[i].size();
//       const char *data = &(*field[i].begin());
//       m_datastore_ptr->getRoot()->createView(view_name + std::to_string(view_count), dtype, 
//                                                          num_elements, (void*)data);
//       view_count++;
//    }
//    return m_datastore_ptr->getRoot()->getView(view_name + "0");
// }

//------------------------------------------------------------------------------
template<typename Dimension>
inline
axom::sidre::View *SidreDataCollection::alloc_view(const std::string &view_name, 
                                                   const Spheral::Field<Dimension, Dim<1>::Vector> &field)
{
   axom::sidre::DataTypeId dtype = field.getAxomType();
   axom::IndexType num_elements = 1;
   int view_count = 0;
   for (u_int i = 0; i < field.size(); i++)
   {
      auto *data = &(*field[i].begin());
      m_datastore_ptr->getRoot()->createView(view_name + std::to_string(view_count), dtype, 
                                                         num_elements, (void*)data);
      view_count++;
   }
   return m_datastore_ptr->getRoot()->getView(view_name + "0");
}

template<typename Dimension>
inline
axom::sidre::View *SidreDataCollection::alloc_view(const std::string &view_name, 
                                                   const Spheral::Field<Dimension, Dim<1>::Vector3d> &field)
{
   axom::sidre::DataTypeId dtype = field.getAxomType();
   axom::IndexType num_elements = 3;
   int view_count = 0;
   for (u_int i = 0; i < field.size(); i++)
   {
      auto *data = &(*field[i].begin());
      m_datastore_ptr->getRoot()->createView(view_name + std::to_string(view_count), dtype, 
                                                         num_elements, (void*)data);
      view_count++;
   }
   return m_datastore_ptr->getRoot()->getView(view_name + "0");
}

template<typename Dimension>
inline
axom::sidre::View *SidreDataCollection::alloc_view(const std::string &view_name, 
                                                   const Spheral::Field<Dimension, Dim<1>::Tensor> &field)
{
   axom::sidre::DataTypeId dtype = field.getAxomType();
   axom::IndexType num_elements = 1;
   int view_count = 0;
   for (u_int i = 0; i < field.size(); i++)
   {
      auto *data = &(*field[i].begin());
      m_datastore_ptr->getRoot()->createView(view_name + std::to_string(view_count), dtype, 
                                                         num_elements, (void*)data);
      view_count++;
   }
   return m_datastore_ptr->getRoot()->getView(view_name + "0");
}

template<typename Dimension>
inline
axom::sidre::View *SidreDataCollection::alloc_view(const std::string &view_name, 
                                                   const Spheral::Field<Dimension, Dim<1>::SymTensor> &field)
{
   axom::sidre::DataTypeId dtype = field.getAxomType();
   axom::IndexType num_elements = 1;
   int view_count = 0;
   for (u_int i = 0; i < field.size(); i++)
   {
      auto *data = &(*field[i].begin());
      m_datastore_ptr->getRoot()->createView(view_name + std::to_string(view_count), dtype, 
                                                         num_elements, (void*)data);
      view_count++;
   }
   return m_datastore_ptr->getRoot()->getView(view_name + "0");
}

template<typename Dimension>
inline
axom::sidre::View *SidreDataCollection::alloc_view(const std::string &view_name, 
                                                   const Spheral::Field<Dimension, Dim<1>::ThirdRankTensor> &field)
{
   axom::sidre::DataTypeId dtype = field.getAxomType();
   axom::IndexType num_elements = 1;
   int view_count = 0;
   for (u_int i = 0; i < field.size(); i++)
   {
      auto *data = &(*field[i].begin());
      m_datastore_ptr->getRoot()->createView(view_name + std::to_string(view_count), dtype, 
                                                         num_elements, (void*)data);
      view_count++;
   }
   return m_datastore_ptr->getRoot()->getView(view_name + "0");
}

template<typename Dimension>
inline
axom::sidre::View *SidreDataCollection::alloc_view(const std::string &view_name, 
                                                   const Spheral::Field<Dimension, Dim<1>::FourthRankTensor> &field)
{
   axom::sidre::DataTypeId dtype = field.getAxomType();
   axom::IndexType num_elements = 1;
   int view_count = 0;
   for (u_int i = 0; i < field.size(); i++)
   {
      auto *data = &(*field[i].begin());
      m_datastore_ptr->getRoot()->createView(view_name + std::to_string(view_count), dtype, 
                                                         num_elements, (void*)data);
      view_count++;
   }
   return m_datastore_ptr->getRoot()->getView(view_name + "0");
}

template<typename Dimension>
inline
axom::sidre::View *SidreDataCollection::alloc_view(const std::string &view_name, 
                                                   const Spheral::Field<Dimension, Dim<1>::FifthRankTensor> &field)
{
   axom::sidre::DataTypeId dtype = field.getAxomType();
   axom::IndexType num_elements = 1;
   int view_count = 0;
   for (u_int i = 0; i < field.size(); i++)
   {
      auto *data = &(*field[i].begin());
      m_datastore_ptr->getRoot()->createView(view_name + std::to_string(view_count), dtype, 
                                                         num_elements, (void*)data);
      view_count++;
   }
   return m_datastore_ptr->getRoot()->getView(view_name + "0");
}

// template<typename Dimension>
// inline
// axom::sidre::View *SidreDataCollection::alloc_view(const std::string &view_name, 
//                                                    const Spheral::Field<Dimension, Dim<1>::FacetedVolume> &field)
// {
//    axom::sidre::DataTypeId dtype = field.getAxomType();
//    axom::IndexType num_elements = 1;
//    int view_count = 0;
//    for (u_int i = 0; i < field.size(); i++)
//    {
//       std::vector<Spheral::Dim<1>::Vector> data = field[i].vertices();
//       std::cout << data.size() << std::endl;
//       for (u_int j = 0; j < data.size(); j++)
//       {
//          auto *data2 = &(*data[i].begin());
//          m_datastore_ptr->getRoot()->createView(view_name + std::to_string(view_count), dtype, 
//                                                          num_elements, (void*)data2);
//       }
      
//       view_count++;
//    }
//    return m_datastore_ptr->getRoot()->getView(view_name + "0");
// }

//------------------------------------------------------------------------------
template<typename Dimension>
inline
axom::sidre::View *SidreDataCollection::alloc_view(const std::string &view_name, 
                                                   const Spheral::Field<Dimension, Dim<2>::Vector> &field)
{
   axom::sidre::DataTypeId dtype = field.getAxomType();
   axom::IndexType num_elements = 2;
   int view_count = 0;
   for (u_int i = 0; i < field.size(); i++)
   {
      auto *data = &(*field[i].begin());
      m_datastore_ptr->getRoot()->createView(view_name + std::to_string(view_count), dtype, 
                                                         num_elements, (void*)data);
      view_count++;
   }
   return m_datastore_ptr->getRoot()->getView(view_name + "0");
}

template<typename Dimension>
inline
axom::sidre::View *SidreDataCollection::alloc_view(const std::string &view_name, 
                                                   const Spheral::Field<Dimension, Dim<2>::Tensor> &field)
{
   axom::sidre::DataTypeId dtype = field.getAxomType();
   axom::IndexType num_elements = 4;
   int view_count = 0;
   for (u_int i = 0; i < field.size(); i++)
   {
      auto *data = &(*field[i].begin());
      m_datastore_ptr->getRoot()->createView(view_name + std::to_string(view_count), dtype, 
                                                         num_elements, (void*)data);
      view_count++;
   }
   return m_datastore_ptr->getRoot()->getView(view_name + "0");
}

template<typename Dimension>
inline
axom::sidre::View *SidreDataCollection::alloc_view(const std::string &view_name, 
                                                   const Spheral::Field<Dimension, Dim<2>::SymTensor> &field)
{
   axom::sidre::DataTypeId dtype = field.getAxomType();
   axom::IndexType num_elements = 2;
   int view_count = 0;
   for (u_int i = 0; i < field.size(); i++)
   {
      auto *data = &(*field[i].begin());
      m_datastore_ptr->getRoot()->createView(view_name + std::to_string(view_count), dtype, 
                                                         num_elements, (void*)data);
      view_count++;
   }
   return m_datastore_ptr->getRoot()->getView(view_name + "0");
}

template<typename Dimension>
inline
axom::sidre::View *SidreDataCollection::alloc_view(const std::string &view_name, 
                                                   const Spheral::Field<Dimension, Dim<2>::ThirdRankTensor> &field)
{
   axom::sidre::DataTypeId dtype = field.getAxomType();
   axom::IndexType num_elements = 1;
   int view_count = 0;
   for (u_int i = 0; i < field.size(); i++)
   {
      auto *data = &(*field[i].begin());
      m_datastore_ptr->getRoot()->createView(view_name + std::to_string(view_count), dtype, 
                                                         num_elements, (void*)data);
      view_count++;
   }
   return m_datastore_ptr->getRoot()->getView(view_name + "0");
}

template<typename Dimension>
inline
axom::sidre::View *SidreDataCollection::alloc_view(const std::string &view_name, 
                                                   const Spheral::Field<Dimension, Dim<2>::FourthRankTensor> &field)
{
   axom::sidre::DataTypeId dtype = field.getAxomType();
   axom::IndexType num_elements = 1;
   int view_count = 0;
   for (u_int i = 0; i < field.size(); i++)
   {
      auto *data = &(*field[i].begin());
      m_datastore_ptr->getRoot()->createView(view_name + std::to_string(view_count), dtype, 
                                                         num_elements, (void*)data);
      view_count++;
   }
   return m_datastore_ptr->getRoot()->getView(view_name + "0");
}

template<typename Dimension>
inline
axom::sidre::View *SidreDataCollection::alloc_view(const std::string &view_name, 
                                                   const Spheral::Field<Dimension, Dim<2>::FifthRankTensor> &field)
{
   axom::sidre::DataTypeId dtype = field.getAxomType();
   axom::IndexType num_elements = 1;
   int view_count = 0;
   for (u_int i = 0; i < field.size(); i++)
   {
      auto *data = &(*field[i].begin());
      m_datastore_ptr->getRoot()->createView(view_name + std::to_string(view_count), dtype, 
                                                         num_elements, (void*)data);
      view_count++;
   }
   return m_datastore_ptr->getRoot()->getView(view_name + "0");
}

//------------------------------------------------------------------------------
template<typename Dimension>
inline
axom::sidre::View *SidreDataCollection::alloc_view(const std::string &view_name, 
                                                   const Spheral::Field<Dimension, Dim<3>::Vector> &field)
{
   axom::sidre::DataTypeId dtype = field.getAxomType();
   axom::IndexType num_elements = 3;
   int view_count = 0;
   for (u_int i = 0; i < field.size(); i++)
   {
      auto *data = &(*field[i].begin());
      m_datastore_ptr->getRoot()->createView(view_name + std::to_string(view_count), dtype, 
                                                         num_elements, (void*)data);
      view_count++;
   }
   return m_datastore_ptr->getRoot()->getView(view_name + "0");
}

template<typename Dimension>
inline
axom::sidre::View *SidreDataCollection::alloc_view(const std::string &view_name, 
                                                   const Spheral::Field<Dimension, Dim<3>::Tensor> &field)
{
   axom::sidre::DataTypeId dtype = field.getAxomType();
   axom::IndexType num_elements = 9;
   int view_count = 0;
   for (u_int i = 0; i < field.size(); i++)
   {
      auto *data = &(*field[i].begin());
      m_datastore_ptr->getRoot()->createView(view_name + std::to_string(view_count), dtype, 
                                                         num_elements, (void*)data);
      view_count++;
   }
   return m_datastore_ptr->getRoot()->getView(view_name + "0");
}

template<typename Dimension>
inline
axom::sidre::View *SidreDataCollection::alloc_view(const std::string &view_name, 
                                                   const Spheral::Field<Dimension, Dim<3>::SymTensor> &field)
{
   axom::sidre::DataTypeId dtype = field.getAxomType();
   axom::IndexType num_elements = 3;
   int view_count = 0;
   for (u_int i = 0; i < field.size(); i++)
   {
      auto *data = &(*field[i].begin());
      m_datastore_ptr->getRoot()->createView(view_name + std::to_string(view_count), dtype, 
                                                         num_elements, (void*)data);
      view_count++;
   }
   return m_datastore_ptr->getRoot()->getView(view_name + "0");
}

template<typename Dimension>
inline
axom::sidre::View *SidreDataCollection::alloc_view(const std::string &view_name, 
                                                   const Spheral::Field<Dimension, Dim<3>::ThirdRankTensor> &field)
{
   axom::sidre::DataTypeId dtype = field.getAxomType();
   axom::IndexType num_elements = 1;
   int view_count = 0;
   for (u_int i = 0; i < field.size(); i++)
   {
      auto *data = &(*field[i].begin());
      m_datastore_ptr->getRoot()->createView(view_name + std::to_string(view_count), dtype, 
                                                         num_elements, (void*)data);
      view_count++;
   }
   return m_datastore_ptr->getRoot()->getView(view_name + "0");
}

template<typename Dimension>
inline
axom::sidre::View *SidreDataCollection::alloc_view(const std::string &view_name, 
                                                   const Spheral::Field<Dimension, Dim<3>::FourthRankTensor> &field)
{
   axom::sidre::DataTypeId dtype = field.getAxomType();
   axom::IndexType num_elements = 1;
   int view_count = 0;
   for (u_int i = 0; i < field.size(); i++)
   {
      auto *data = &(*field[i].begin());
      m_datastore_ptr->getRoot()->createView(view_name + std::to_string(view_count), dtype, 
                                                         num_elements, (void*)data);
      view_count++;
   }
   return m_datastore_ptr->getRoot()->getView(view_name + "0");
}

template<typename Dimension>
inline
axom::sidre::View *SidreDataCollection::alloc_view(const std::string &view_name, 
                                                   const Spheral::Field<Dimension, Dim<3>::FifthRankTensor> &field)
{
   axom::sidre::DataTypeId dtype = field.getAxomType();
   axom::IndexType num_elements = 1;
   int view_count = 0;
   for (u_int i = 0; i < field.size(); i++)
   {
      auto *data = &(*field[i].begin());
      m_datastore_ptr->getRoot()->createView(view_name + std::to_string(view_count), dtype, 
                                                         num_elements, (void*)data);
      view_count++;
   }
   return m_datastore_ptr->getRoot()->getView(view_name + "0");
}


}
