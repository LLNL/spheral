#ifndef SidreDataCollection_HH
#define SidreDataCollection_HH

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
    axom::sidre::View *alloc_view(const std::string &view_name, 
                                  const Spheral::Field<Dimension, DataType> &field);
    
    void printDataStore() {m_datastore_ptr->getRoot()->print();};
private:
    axom::sidre::DataStore *m_datastore_ptr;
};

}

#include "Utilities/SidreDataCollectionInline.hh"

#endif
