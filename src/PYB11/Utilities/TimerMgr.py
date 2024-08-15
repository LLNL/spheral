#-------------------------------------------------------------------------------
# TimerMgr class
#-------------------------------------------------------------------------------
from PYB11Generator import *

@PYB11singleton
class TimerMgr:

    "Singleton wrapper for CaliperManager. Access through TimerMgr.instance(), ie TimerMgr.instance().start()."

    @PYB11static
    def instance(self):
        "Access the singleton instance of the timer manager"
        return "TimerMgr&"

    def timer_start(self, region_name = "std::string"):
        "Start custom region Caliper timer, must have corresponding timer_end call"
        return "void"

    def timer_end(self, region_name = "std::string"):
        "End custom region Caliper timer"
        return "void"

    def add(self, config_str = "std::string"):
        "Add a Caliper configuration"
        return "void"

    def default_start(self, testname = "std::string"):
        "Set the spot Caliper configuration and start the manager"
        return "void"

    def start(self):
        "Start the Caliper configuration manager"
        return "void"

    def stop(self):
        "Stop the Caliper configuration manager"
        return "void"

    def fini(self):
        "Flush the Caliper configuration manager"
        return "void"
