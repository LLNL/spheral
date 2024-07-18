#-------------------------------------------------------------------------------
# Adiak utilities
#-------------------------------------------------------------------------------
from PYB11Generator import *

# This is defined in the Utilities_PYB11.py preamble
@PYB11cppname("Spheral::spheral_adiak_init")
def adiak_init():
    "Initialize Adiak and run collect_all"
    return "void"

@PYB11cppname("adiak::fini")
def adiak_fini():
    "Finalize Adiak"
    return "void"

@PYB11cppname("adiak::collect_all")
def adiak_collect_all():
    "Collect all default Adiak metadata"
    return "void"

adiak_categories = PYB11enum(("unset", "all", "general", "performance", "control"),
                             doc="Enum of Adiak categories")

@PYB11cppname("adiak::value")
@PYB11template("ValueType")
def adiak_value(name = "std::string",
                value = "%(ValueType)s",
                category = ("int",
                            "adiak_categories::general"),
                subcategory = ("std::string", '""')):
    "Set a single value in Adiak with a given name"
    return "bool"

@PYB11cppname("adiak::value")
@PYB11pyname("adiak_value")
@PYB11template("ValueType")
def adiak_value2(name = "std::string",
                 value = "%(ValueType)s",
                 value2 = "%(ValueType)s",
                 category = ("int",
                             "adiak_categories::general"),
                 subcategory = ("std::string", '""')):
    "Set a pair of values in Adiak with a given name"
    return "bool"
