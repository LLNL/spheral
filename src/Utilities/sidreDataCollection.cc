#include "sidreDataCollection.hh"


namespace sidre = axom::sidre;


namespace Spheral {

SidreDataCollection::SidreDataCollection()
{
   m_datastore_ptr = new sidre::DataStore();
   sidre::Group *root = m_datastore_ptr->getRoot();
}

SidreDataCollection::~SidreDataCollection()
{
   delete m_datastore_ptr;
}

}
