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
   auto *data = &(*field.begin());
   axom::sidre::View *v = m_datastore_ptr->getRoot()->createView(view_name, dtype,
                                                        num_elements, (void*)data);
   return v;
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
      auto *data = &(*field[i].begin());
      m_datastore_ptr->getRoot()->createView(view_name + std::to_string(view_count), dtype, 
                                                         num_elements, (void*)data);
      view_count++;
   }
   return m_datastore_ptr->getRoot()->getView(view_name + "0");
}

// //------------------------------------------------------------------------------
// template<typename Dimension, typename DataType>
//     axom::sidre::View *alloc_view(const std::string &view_name,
//                                   const Spheral::Field<Dimension, std::tuple<DataType, DataType, DataType>> &field)
// {
//    axom::sidre::DataTypeId dtype = field.getAxomType();
//    axom::IndexType num_elements = field[i].size();
//    int view_count = 0;
//    for (u_int i = 0; i < field.size(); i++)
//    {
//       auto *data = &(*field[i].begin());
//       m_datastore_ptr->getRoot()->createView(view_name + std::to_string(view_count), dtype, 
//                                                          num_elements, (void*)data);
//       view_count++;
//    }
//    return m_datastore_ptr->getRoot()->getView(view_name + "0");
// }

// //------------------------------------------------------------------------------
// template<typename Dimension, typename DataType>
//     axom::sidre::View *alloc_view(const std::string &view_name,
//                                   const Spheral::Field<Dimension, std::tuple<DataType, DataType, DataType, DataType>> &field)
// {
//    axom::sidre::DataTypeId dtype = field.getAxomType();
//    axom::IndexType num_elements = field[i].size();
//    int view_count = 0;
//    for (u_int i = 0; i < field.size(); i++)
//    {
//       auto *data = &(*field[i].begin());
//       m_datastore_ptr->getRoot()->createView(view_name + std::to_string(view_count), dtype, 
//                                                          num_elements, (void*)data);
//       view_count++;
//    }
//    return m_datastore_ptr->getRoot()->getView(view_name + "0");
// }

// //------------------------------------------------------------------------------
// template<typename Dimension, typename DataType>
//     axom::sidre::View *alloc_view(const std::string &view_name,
//                                   const Spheral::Field<Dimension, std::tuple<DataType, DataType, DataType, DataType, DataType>> &field)
// {
//    axom::sidre::DataTypeId dtype = field.getAxomType();
//    axom::IndexType num_elements = field[i].size();
//    int view_count = 0;
//    for (u_int i = 0; i < field.size(); i++)
//    {
//       auto *data = &(*field[i].begin());
//       m_datastore_ptr->getRoot()->createView(view_name + std::to_string(view_count), dtype, 
//                                                          num_elements, (void*)data);
//       view_count++;
//    }
//    return m_datastore_ptr->getRoot()->getView(view_name + "0");
}



}
