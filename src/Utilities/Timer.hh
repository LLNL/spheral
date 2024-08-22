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
  bool started = false;
  std::string caliperFilename = "";
  std::string caliperConfig = "";
public:
  static TimerMgr& instance() {
    static TimerMgr theInstance;
    return theInstance;
  }
  static void timer_start(std::string regionName) {
    TIME_BEGIN(regionName.c_str());
  }
  static void timer_end(std::string regionName) {
    TIME_END(regionName.c_str());
  }
  static bool is_started() {
    return instance().started;
  }
  static std::string get_config() {
    return instance().caliperConfig;
  }
  static std::string get_filename() {
    return instance().caliperFilename;
  }
#ifdef TIMER
private:
  cali::ConfigManager cali_mgr;
public:
  static void add(std::string config_str) {
    bool test = instance().cali_mgr.add(config_str.c_str());
    VERIFY2(test, instance().cali_mgr.error_msg());
    instance().caliperConfig += config_str;
  }
  static void default_start(std::string testname) {
    if (!testname.empty()) {
      std::string default_config = "spot,mem.highwatermark,output=" + testname + ".cali";
      instance().caliperFilename = testname + ".cali";
      add(default_config);
      start();
    } else if (Spheral::Process::getRank() == 0) {
      std::cout << "WARNING: Caliper test name is empty, "
                << "no Caliper configuration started" << std::endl;
    }
  }
  static void start() {
    instance().cali_mgr.start();
    instance().started = true;
  }
  static void stop() {
    instance().cali_mgr.stop();
  }
  static void fini() {
    instance().cali_mgr.flush();
    instance().started = false;
  }
#else
  static void default_start(std::string) {
  }
  static void add(std::string) {
  }
  static void start() {
  }
  static void stop() {
  }
  static void fini() {
  }
#endif
};
}
#endif
