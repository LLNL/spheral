#include "testSidreStorage.hh"

TYPED_TEST(SidreDataCollectionTestNew, SidreAllocViewTest)
{
  Spheral::SidreDataCollection myData;
  size_t n = 10;

  SpheralTestNodeList<TypeParam> makeNodeList("test bed", n, 0);
  SpheralTestField<TypeParam> testField("test field", makeNodeList);

  initField<TypeParam>(testField, n);

  axom::sidre::Group *myTemp = myData.sidreStoreField("SidreTest", testField);

  testSidreData<TypeParam>(testField, myTemp);
}
