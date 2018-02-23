//------------------------------------------------------------------------------
// Wrapper methods to control the EasyProfiler timers.
// https://github.com/yse/easy_profiler
//------------------------------------------------------------------------------
#ifndef __Spheral_EasyProfileer_helpers__
#define __Spheral_EasyProfileer_helpers__

#include "Utilities/Process.hh"
#include "easy/profiler.h"
#include <string>

//------------------------------------------------------------------------------
// EasyProfilerStart
//------------------------------------------------------------------------------
inline
void EasyProfilerStart() {
  EASY_PROFILER_ENABLE;
}

//------------------------------------------------------------------------------
// EasyProfilerDump
//------------------------------------------------------------------------------
inline
void EasyProfilerDump(const std::string basename) {
  const auto rank = Process::getRank();
  const std::string filename = basename + "_" + std::to_string(rank) + ".easyprof";
  profiler::dumpBlocksToFile(filename.c_str());
}

#endif
