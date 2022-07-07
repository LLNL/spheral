#include <caliper/cali.h>
#include <iostream>

int main()
{
  // CALI_CXX_MARK_FUNCTION;
  CALI_MARK_BEGIN("area");
  std::cout << "Hello World!" << std::endl;
  CALI_MARK_END("area");

  int a = 12398745;

  CALI_CXX_MARK_LOOP_BEGIN(mainloop_id, "mainloop");
  for (int i = 0; i < 100; ++i)
  {
    CALI_CXX_MARK_LOOP_ITERATION(mainloop_id, i);
    a += 13;
  }
  CALI_CXX_MARK_LOOP_END(mainloop_id);
  return 0;
}