//---------------------------------Spheral++----------------------------------//
// Timer macros and Caliper ConfigManger singleton
//
//----------------------------------------------------------------------------//
#ifndef __Spheral_Timer__
#define __Spheral_Timer__
// If TIMER is defined then we want timer functionality
#ifdef TIMER
#include "caliper/cali.h"
#include "caliper/cali-manager.h"

#define TIME_FUNCTION CALI_CXX_MARK_FUNCTION
#define TIME_BEGIN(regionName) CALI_MARK_BEGIN(regionName)
#define TIME_END(regionName) CALI_MARK_END(regionName)

#else // TIMER
// Stub TIME macros, when TIMER is off

#define TIME_FUNCTION
#define TIME_BEGIN(regionName)
#define TIME_END(regionName)

#endif // TIMER
// Note: This class is initialized in
// SimulationControl/SpheralOptionParser.py
namespace Spheral {
class TimerMgr {
private:
  TimerMgr() = default;
  ~TimerMgr() { }
  TimerMgr(const TimerMgr&) = delete;
  TimerMgr& operator=(const TimerMgr&) = delete;
public:
  static TimerMgr& instance() {
    static TimerMgr theInstance;
    return theInstance;
  }
  void timer_start(std::string regionName) {
    TIME_BEGIN(regionName.c_str());
  }
  void timer_end(std::string regionName) {
    TIME_END(regionName.c_str());
  }
#ifdef TIMER
private:
  cali::ConfigManager cali_mgr;
public:
  void add(std::string config_str) {
    bool test = cali_mgr.add(config_str.c_str());
    VERIFY2(test, cali_mgr.error_msg());
  }
  void default_start(std::string testname) {
    if (!testname.empty()) {
      std::string default_config = "spot,output=" + testname + ".cali,mem.highwatermark";
      add(default_config);
      start();
    } else if (Spheral::Process::getRank() == 0) {
      std::cout << "WARNING: Caliper test name is empty, "
                << "no Caliper configuration started" << std::endl;
    }
  }
  void start() {
    cali_mgr.start();
  }
  void stop() {
    cali_mgr.stop();
  }
  void fini() {
    cali_mgr.flush();
  }
#else
  void default_start(std::string) {
  }
  void add(std::string) {
  }
  void start() {
  }
  void stop() {
  }
  void fini() {
  }
#endif
};
}
#endif
