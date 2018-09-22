"""
Spheral CXXTypes module.

Provides access to fundamental C++ types from python that are not Spheral specific.
"""

from PYB11Generator import *

# preamble = """
# typedef std::pair<double, double> pair_double_double;
# typedef std::pair<double, std::string> pair_double_string;
# typedef std::pair<unsigned, unsigned> pair_unsigned_unsigned;
# typedef std::pair<uint64_t, uint64_t> pair_ULL_ULL;
# typedef std::pair<std::string, std::string> pair_string_string;
# """

# Include files.
includes = ["<vector>",
            "<map>",
            "<set>",
            "<string>"]

# std::vector
vector_of_char     = PYB11_bind_vector("char", opaque=True)
vector_of_unsigned = PYB11_bind_vector("unsigned", opaque=True)
vector_of_ULL      = PYB11_bind_vector("uint64_t", opaque=True)
vector_of_int      = PYB11_bind_vector("int", opaque=True)
vector_of_float    = PYB11_bind_vector("float", opaque=True)
vector_of_double   = PYB11_bind_vector("double", opaque=True)
vector_of_string   = PYB11_bind_vector("std::string", opaque=True)

# std::vector<std::vector>
vector_of_vector_of_char     = PYB11_bind_vector("std::vector<char>", opaque=True)
vector_of_vector_of_unsigned = PYB11_bind_vector("std::vector<unsigned>", opaque=True)
vector_of_vector_of_ULL      = PYB11_bind_vector("std::vector<uint64_t>", opaque=True)
vector_of_vector_of_int      = PYB11_bind_vector("std::vector<int>", opaque=True)
vector_of_vector_of_float    = PYB11_bind_vector("std::vector<float>", opaque=True)
vector_of_vector_of_double   = PYB11_bind_vector("std::vector<double>", opaque=True)
vector_of_vector_of_string   = PYB11_bind_vector("std::vector<std::string>", opaque=True)

# std::vector<pair<>>
vector_of_pair_double_double     = PYB11_bind_vector("pair<double,_double>", opaque=True)
# vector_of_pair_double_string     = PYB11_bind_vector("pair_double_string", opaque=True)
# vector_of_pair_unsigned_unsigned = PYB11_bind_vector("pair_unsigned_unsigned", opaque=True)
# vector_of_pair_ULL_ULL           = PYB11_bind_vector("pair_ULL_ULL", opaque=True)
# vector_of_pair_string_string     = PYB11_bind_vector("pair_string_string", opaque=True)

# std::map
map_string_double = PYB11_bind_map("std::string", "double", opaque=True)
map_int_string    = PYB11_bind_map("int", "std::string", opaque=True)
