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

}
