#include "axom/sidre.hpp"
#include "Field/Field.hh"

namespace Spheral
{

template<typename Dimension, typename DataType>
inline
/* sidre::View * */ void SidreDataCollection::alloc_view(const std::string &view_name, 
                                             const Spheral::Field<Dimension, DataType> &field)
{
   axom::sidre::DataTypeId dtype = field.getAxomType();
   axom::IndexType num_elements = field.numElements();
   void *data = field.allValues().data();
   axom::sidre::View *v = m_datastore_ptr->getRoot()->createView(view_name, dtype,
                                                           num_elements, data);
}

}
