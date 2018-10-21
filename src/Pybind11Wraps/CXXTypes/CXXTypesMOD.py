"""
Spheral CXXTypes module.

Provides access to fundamental C++ types from python that are not Spheral specific.
"""

from PYB11Generator import *

# Include files.
includes = ["<vector>",
            "<map>",
            "<set>",
            "<string>"]

# std::vector
vector_of_char     = PYB11_bind_vector("char", opaque=True, local=False)
vector_of_unsigned = PYB11_bind_vector("unsigned", opaque=True, local=False)
vector_of_ULL      = PYB11_bind_vector("uint64_t", opaque=True, local=False)
vector_of_int      = PYB11_bind_vector("int", opaque=True, local=False)
vector_of_float    = PYB11_bind_vector("float", opaque=True, local=False)
vector_of_double   = PYB11_bind_vector("double", opaque=True, local=False)
vector_of_string   = PYB11_bind_vector("std::string", opaque=True, local=False)

# std::vector<std::vector>
vector_of_vector_of_char     = PYB11_bind_vector("std::vector<char>", opaque=True, local=False)
vector_of_vector_of_unsigned = PYB11_bind_vector("std::vector<unsigned>", opaque=True, local=False)
vector_of_vector_of_ULL      = PYB11_bind_vector("std::vector<uint64_t>", opaque=True, local=False)
vector_of_vector_of_int      = PYB11_bind_vector("std::vector<int>", opaque=True, local=False)
vector_of_vector_of_float    = PYB11_bind_vector("std::vector<float>", opaque=True, local=False)
vector_of_vector_of_double   = PYB11_bind_vector("std::vector<double>", opaque=True, local=False)
vector_of_vector_of_string   = PYB11_bind_vector("std::vector<std::string>", opaque=True, local=False)

# std::vector<pair<>>
vector_of_pair_double_double     = PYB11_bind_vector("std::pair<double, double>", opaque=True, local=False)
vector_of_pair_double_string     = PYB11_bind_vector("std::pair<double, std::string>", opaque=True, local=False)
vector_of_pair_unsigned_unsigned = PYB11_bind_vector("std::pair<unsigned, unsigned>", opaque=True, local=False)
vector_of_pair_ULL_ULL           = PYB11_bind_vector("std::pair<uint64_t, uint64_t>", opaque=True, local=False)
vector_of_pair_string_string     = PYB11_bind_vector("std::pair<std::string, std::string>", opaque=True, local=False)

# std::map
map_string_double = PYB11_bind_map("std::string", "double", opaque=True, local=False)
map_int_string    = PYB11_bind_map("int", "std::string", opaque=True, local=False)
