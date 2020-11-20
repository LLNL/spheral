#include "gtest/gtest.h"

#include "Field/Field.hh"
#include "Utilities/SidreDataCollection.hh"
#include "Geometry/Dimension.hh"

namespace Spheral
{
    TEST(basic, example)
    {
        Field<Dim<1>, int> testField ("test field");
        //ASSERT_EQ(2,2);
    }


}

