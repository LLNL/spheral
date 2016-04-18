//------------------------------------------------------------------------------
// timingUtilities
//
// A set of inline helper methods to encapsulate how we do timing.
//
// JMO:  Tue Dec  9 10:31:14 PST 2008
//------------------------------------------------------------------------------
#define BOOST_DATE_TIME_POSIX_TIME_STD_CONFIG     // Enable nanosecond timings.
#include <boost/date_time/posix_time/posix_time.hpp>

namespace Spheral {
//------------------------------------------------------------------------------
// Get the current clock time.
//------------------------------------------------------------------------------
struct Timing {
  typedef boost::posix_time::ptime Time;
  typedef boost::posix_time::time_duration duration;
  static Time currentTime() { return boost::posix_time::microsec_clock::local_time(); }
  static double convertToSeconds(const duration& delta) { return double(delta.total_nanoseconds())/1e9; }
  static double difference(const Time& t1, const Time& t2) { return convertToSeconds(t2 - t1); }
};
}
