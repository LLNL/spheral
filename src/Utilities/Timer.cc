#include "Utilities/Timer.hh"
#ifdef TIMER
#include "Utilities/BuildData.hh"
#include "Utilities/OpenMP_wrapper.hh"
#include "adiak.hpp"
#include <omp.h>
#endif // TIMER

#include <string>

namespace Spheral {

void initAdiakData(std::string test_name,
                   std::string spheral_branch,
                   std::string spheral_short_commit,
                   int problem_size)
{
#ifdef TIMER
  adiak::init(NULL);
  adiak::launchdate();
  adiak::executablepath();
  adiak::cmdline();
  adiak::clustername();
  adiak::jobsize();
  adiak::value("test_name", test_name);
  adiak::value("problem_size", problem_size);
  adiak::value("spheral_branch", spheral_branch);
  adiak::value("spheral_short_commit", spheral_short_commit);
  adiak::value("omp_num_threads", omp_get_max_threads());
  adiak::value("cxx_compiler_id", BuildData::cxx_compiler_id);
  adiak::fini();
#endif // TIMER
}

} //  namespace Spheral
