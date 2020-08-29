#pragma once

#include "axom/sidre.hpp"
#include "Field/Field.hh"

class SidreDataCollection
{
public:
    SidreDataCollection();
    ~SidreDataCollection();

    template<typename Dimension, typename DataType>
    axom::sidre::View *alloc_view(const std::string &view_name, 
                                  const Spheral::Field<Dimension, DataType> &field);
private:
    axom::sidre::DataStore *m_datastore_ptr;
};
