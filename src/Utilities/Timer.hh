// If TIMER is defined then we want timer functionality
#ifdef TIMER
#include "caliper/cali.h"

#define TIME_FUNCTION CALI_CXX_MARK_FUNCTION
#define TIME_BEGIN(regionName) CALI_MARK_BEGIN(regionName)
#define TIME_END(regionName) CALI_MARK_END(regionName)

#else // TIMER
// Stub TIME macros, when TIME is zero
#define TIME_FUNCTION
#define TIME_BEGIN(regionName)
#define TIME_END(regionName)

#endif // TIMER