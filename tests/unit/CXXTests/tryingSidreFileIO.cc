#include <iostream>
#include <vector>
#include <memory>
#include <string>

#include "axom/sidre.hpp"
#include "axom/slic.hpp"

using namespace std;


int main()
{
  axom::sidre::DataStore mFilePtr;
  axom::sidre::DataStore ds;
  axom::sidre::Group* root2 = mFilePtr.getRoot();
  axom::sidre::Group* root = ds.getRoot();
  string pathName = "path";
  string fileName = "file";

  int size = 20;
  vector<int> data;
  vector<int> value;

  for (int i = 0; i < size; ++i)
  {
    data.push_back((i + 1) * 10);
  }

  root2->createView(pathName, axom::sidre::INT_ID, data.size(), (void*)(&(*data.begin())));
  root2->print();
  std::cout << std::endl;

  root2->save(fileName, "sidre_hdf5");

  root->load(fileName);

  root->print();

  SLIC_ASSERT(root->hasChildView(pathName));
  axom::sidre::View* view = root->getView(pathName);
  SLIC_ASSERT(view->isExternal());

  int valSize = view->getNumElements();
  value.resize(valSize);
  std::cout << value.size() << "\n";

  view->setExternalDataPtr(static_cast<void*>(&value[0]));
  root->loadExternalData(fileName);

  root->print();


  return 0;
}