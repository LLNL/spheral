#ifndef sidreDataCollection_HH
#define sidreDataCollection_HH

#include "axom/sidre.hpp"
#include "Field/Field.hh"

namespace Spheral
{

class SidreDataCollection
{
public:
    SidreDataCollection();
    ~SidreDataCollection();

    template<typename Dimension, typename DataType>
    /* axom::sidre::View * */ void alloc_view(const std::string &view_name, 
                                  const Spheral::Field<Dimension, DataType> &field);
private:
    axom::sidre::DataStore *m_datastore_ptr;
};

}

#include "Utilities/sidreDataCollectionInline.hh"

#endif
