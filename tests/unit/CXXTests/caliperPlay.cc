#include <caliper/cali.h>
#include <iostream>

int main()
{
  // CALI_CXX_MARK_FUNCTION;
  CALI_MARK_BEGIN("area");
  std::cout << "Hello World!" << std::endl;
  CALI_MARK_END("area");
  return 0;
}