#-------------------------------------------------------------------------------
# Helper to define methods for pickling/unpickling Spheral types.
# This can be important with pyMPI, since pyMPI relies on picking data types
# for communication!
#-------------------------------------------------------------------------------
import copy_reg, pickle
from SpheralModules import *
from SpheralModules.Spheral import *
from SpheralModules.Spheral.NodeSpace import *
from SpheralModules.Spheral.FieldSpace import *
from SpheralModules.Spheral.DataBaseSpace import *
from SpheralModules.Spheral.FileIOSpace import *
from SpheralModules.Spheral.ArtificialViscositySpace import *
from SpheralModules.Spheral.DataOutput import *
from SpheralModules.Spheral.KernelSpace import *
from SpheralModules.Spheral.NeighborSpace import *
from SpheralModules.Spheral.Material import *
from SpheralModules.Spheral.BoundarySpace import *
from SpheralModules.Spheral.PhysicsSpace import *
from SpheralModules.Spheral.GravitySpace import *
from SpheralModules.Spheral.IntegratorSpace import *

## from SpheralModules.Spheral import Vector1d, Vector2d, Vector3d, \
##                                    Tensor1d, Tensor2d, Tensor3d, \
##                                    SymTensor1d, SymTensor2d, SymTensor3d, \
##                                    ThirdRankTensor1d, ThirdRankTensor2d, ThirdRankTensor3d

#-------------------------------------------------------------------------------
# Vector1d
#-------------------------------------------------------------------------------
def construct_Vector1d(x):
    return Vector1d(x)

def reduce_Vector1d(obj):
    return construct_Vector1d, (obj.x,)

copy_reg.pickle(type(Vector1d()), reduce_Vector1d, construct_Vector1d)

#-------------------------------------------------------------------------------
# Vector2d
#-------------------------------------------------------------------------------
def construct_Vector2d(x, y):
    return Vector2d(x, y)

def reduce_Vector2d(obj):
    return construct_Vector2d, (obj.x, obj.y)

copy_reg.pickle(type(Vector2d()), reduce_Vector2d, construct_Vector2d)

#-------------------------------------------------------------------------------
# Vector3d
#-------------------------------------------------------------------------------
def construct_Vector3d(x, y, z):
    return Vector3d(x, y, z)

def reduce_Vector3d(obj):
    return construct_Vector3d, (obj.x, obj.y, obj.z)

copy_reg.pickle(type(Vector3d()), reduce_Vector3d, construct_Vector3d)

#-------------------------------------------------------------------------------
# Tensor1d
#-------------------------------------------------------------------------------
def construct_Tensor1d(xx):
    return Tensor1d(xx)

def reduce_Tensor1d(obj):
    return construct_Tensor1d, (obj.xx,)

copy_reg.pickle(type(Tensor1d()), reduce_Tensor1d, construct_Tensor1d)

#-------------------------------------------------------------------------------
# Tensor2d
#-------------------------------------------------------------------------------
def construct_Tensor2d(xx, xy,
                       yx, yy):
    return Tensor2d(xx, xy,
                    yx, yy)

def reduce_Tensor2d(obj):
    return construct_Tensor2d, (obj.xx, obj.xy,
                                obj.yx, obj.yy)

copy_reg.pickle(type(Tensor2d()), reduce_Tensor2d, construct_Tensor2d)

#-------------------------------------------------------------------------------
# Tensor3d
#-------------------------------------------------------------------------------
def construct_Tensor3d(xx, xy, xz,
                       yx, yy, yz,
                       zx, zy, zz):
    return Tensor3d(xx, xy, xz,
                    yx, yy, yz,
                    zx, zy, zz)

def reduce_Tensor3d(obj):
    return construct_Tensor3d, (obj.xx, obj.xy, obj.xz,
                                obj.yx, obj.yy, obj.yz,
                                obj.zx, obj.zy, obj.zz)

copy_reg.pickle(type(Tensor3d()), reduce_Tensor3d, construct_Tensor3d)

#-------------------------------------------------------------------------------
# SymTensor1d
#-------------------------------------------------------------------------------
def construct_SymTensor1d(xx):
    return SymTensor1d(xx)

def reduce_SymTensor1d(obj):
    return construct_SymTensor1d, (obj.xx,)

copy_reg.pickle(type(SymTensor1d()), reduce_SymTensor1d, construct_SymTensor1d)

#-------------------------------------------------------------------------------
# SymTensor2d
#-------------------------------------------------------------------------------
def construct_SymTensor2d(xx, xy,
                          yx, yy):
    return SymTensor2d(xx, xy,
                       yx, yy)

def reduce_SymTensor2d(obj):
    return construct_SymTensor2d, (obj.xx, obj.xy,
                                   obj.yx, obj.yy)

copy_reg.pickle(type(SymTensor2d()), reduce_SymTensor2d, construct_SymTensor2d)

#-------------------------------------------------------------------------------
# SymTensor3d
#-------------------------------------------------------------------------------
def construct_SymTensor3d(xx, xy, xz,
                          yx, yy, yz,
                          zx, zy, zz):
    return SymTensor3d(xx, xy, xz,
                       yx, yy, yz,
                       zx, zy, zz)

def reduce_SymTensor3d(obj):
    return construct_SymTensor3d, (obj.xx, obj.xy, obj.xz,
                                   obj.yx, obj.yy, obj.yz,
                                   obj.zx, obj.zy, obj.zz)

copy_reg.pickle(type(SymTensor3d()), reduce_SymTensor3d, construct_SymTensor3d)

#-------------------------------------------------------------------------------
# ThirdRankTensor1d
#-------------------------------------------------------------------------------
def construct_ThirdRankTensor1d(x00):
    result = ThirdRankTensor1d()
    for i in xrange(1):
        result[i] = eval("x%02i" % i)
    return result

def reduce_ThirdRankTensor1d(obj):
    return construct_ThirdRankTensor1d, tuple(obj[i] for i in xrange(1))

copy_reg.pickle(type(ThirdRankTensor1d()), reduce_ThirdRankTensor1d, construct_ThirdRankTensor1d)

#-------------------------------------------------------------------------------
# ThirdRankTensor2d
#-------------------------------------------------------------------------------
def construct_ThirdRankTensor2d(x00, x01, x02, x03, x04, x05, x06, x07):
    result = ThirdRankTensor2d()
    for i in xrange(8):
        result[i] = eval("x%02i" % i)
    return result

def reduce_ThirdRankTensor2d(obj):
    return construct_ThirdRankTensor2d, tuple(obj[i] for i in xrange(8))

copy_reg.pickle(type(ThirdRankTensor2d()), reduce_ThirdRankTensor2d, construct_ThirdRankTensor2d)

#-------------------------------------------------------------------------------
# ThirdRankTensor3d
#-------------------------------------------------------------------------------
def construct_ThirdRankTensor3d(x00, x01, x02, x03, x04, x05, x06, x07, x08, x09,
                                x10, x11, x12, x13, x14, x15, x16, x17, x18, x19,
                                x20, x21, x22, x23, x24, x25, x26):
    result = ThirdRankTensor3d()
    for i in xrange(27):
        result[i] = eval("x%02i" % i)
    return result

def reduce_ThirdRankTensor3d(obj):
    return construct_ThirdRankTensor3d, tuple(obj[i] for i in xrange(27))

copy_reg.pickle(type(ThirdRankTensor3d()), reduce_ThirdRankTensor3d, construct_ThirdRankTensor3d)

#-------------------------------------------------------------------------------
# Box1d
#-------------------------------------------------------------------------------
def construct_Box1d(encoded_string):
    return unpackElementBox1d(encoded_string)

def reduce_Box1d(obj):
    return construct_Box1d, (packElementBox1d(obj),)

copy_reg.pickle(type(Box1d()), reduce_Box1d, construct_Box1d)

#-------------------------------------------------------------------------------
# Polygon
#-------------------------------------------------------------------------------
def construct_Polygon(encoded_string):
    return unpackElementPolygon(encoded_string)

def reduce_Polygon(obj):
    return construct_Polygon, (packElementPolygon(obj),)

copy_reg.pickle(type(Polygon()), reduce_Polygon, construct_Polygon)

#-------------------------------------------------------------------------------
# Polyhedron
#-------------------------------------------------------------------------------
def construct_Polyhedron(encoded_string):
    return unpackElementPolyhedron(encoded_string)

def reduce_Polyhedron(obj):
    return construct_Polyhedron, (packElementPolyhedron(obj),)

copy_reg.pickle(type(Polyhedron()), reduce_Polyhedron, construct_Polyhedron)

#------------------------------------------------------------------------------
# std::vectors
#------------------------------------------------------------------------------
# Sadly this nifty bit of exec code doesn't seem to work, so we have to crunch
# out the individual versions explicitly.

## vector_template = """
## def construct_vector_of_%(value_type)s(encoded_string):
##     return string2vector_of_%(value_type)s(encoded_string)[0]
## def reduce_vector_of_%(value_type)s(obj):
##     return construct_vector_of_%(value_type)s, (vector2string(obj, 20),)
## copy_reg.constructor(construct_vector_of_%(value_type)s)
## copy_reg.pickle(vector_of_%(value_type)s, reduce_vector_of_%(value_type)s, construct_vector_of_%(value_type)s)
## """

## vector_value_types = ("int", "double", "string",
##                       "Vector1d", "Tensor1d", "SymTensor1d", "ThirdRankTensor1d",
##                       "Vector2d", "Tensor2d", "SymTensor2d", "ThirdRankTensor2d",
##                       "Vector3d", "Tensor3d", "SymTensor3d", "ThirdRankTensor3d")

def construct_vector_of_int(encoded_string):
    return string2vector_of_int(encoded_string)
def reduce_vector_of_int(obj):
    return construct_vector_of_int, (vector2string(obj, 20),)
copy_reg.constructor(construct_vector_of_int)
copy_reg.pickle(vector_of_int, reduce_vector_of_int, construct_vector_of_int)

def construct_vector_of_ULL(encoded_string):
    return string2vector_of_ULL(encoded_string)
def reduce_vector_of_ULL(obj):
    return construct_vector_of_ULL, (vector2string(obj, 20),)
copy_reg.constructor(construct_vector_of_ULL)
copy_reg.pickle(vector_of_ULL, reduce_vector_of_ULL, construct_vector_of_ULL)

def construct_vector_of_double(encoded_string):
    return string2vector_of_double(encoded_string)
def reduce_vector_of_double(obj):
    return construct_vector_of_double, (vector2string(obj, 20),)
copy_reg.constructor(construct_vector_of_double)
copy_reg.pickle(vector_of_double, reduce_vector_of_double, construct_vector_of_double)

def construct_vector_of_string(encoded_string):
    return string2vector_of_string(encoded_string)
def reduce_vector_of_string(obj):
    return construct_vector_of_string, (vector2string(obj, 20),)
copy_reg.constructor(construct_vector_of_string)
copy_reg.pickle(vector_of_string, reduce_vector_of_string, construct_vector_of_string)

def construct_vector_of_Vector1d(encoded_string):
    return string2vector_of_Vector1d(encoded_string)
def reduce_vector_of_Vector1d(obj):
    return construct_vector_of_Vector1d, (vector2string(obj, 20),)
copy_reg.constructor(construct_vector_of_Vector1d)
copy_reg.pickle(vector_of_Vector1d, reduce_vector_of_Vector1d, construct_vector_of_Vector1d)

def construct_vector_of_Tensor1d(encoded_string):
    return string2vector_of_Tensor1d(encoded_string)
def reduce_vector_of_Tensor1d(obj):
    return construct_vector_of_Tensor1d, (vector2string(obj, 20),)
copy_reg.constructor(construct_vector_of_Tensor1d)
copy_reg.pickle(vector_of_Tensor1d, reduce_vector_of_Tensor1d, construct_vector_of_Tensor1d)

def construct_vector_of_SymTensor1d(encoded_string):
    return string2vector_of_SymTensor1d(encoded_string)
def reduce_vector_of_SymTensor1d(obj):
    return construct_vector_of_SymTensor1d, (vector2string(obj, 20),)
copy_reg.constructor(construct_vector_of_SymTensor1d)
copy_reg.pickle(vector_of_SymTensor1d, reduce_vector_of_SymTensor1d, construct_vector_of_SymTensor1d)

def construct_vector_of_ThirdRankTensor1d(encoded_string):
    return string2vector_of_ThirdRankTensor1d(encoded_string)
def reduce_vector_of_ThirdRankTensor1d(obj):
    return construct_vector_of_ThirdRankTensor1d, (vector2string(obj, 20),)
copy_reg.constructor(construct_vector_of_ThirdRankTensor1d)
copy_reg.pickle(vector_of_ThirdRankTensor1d, reduce_vector_of_ThirdRankTensor1d, construct_vector_of_ThirdRankTensor1d)

def construct_vector_of_Vector2d(encoded_string):
    return string2vector_of_Vector2d(encoded_string)
def reduce_vector_of_Vector2d(obj):
    return construct_vector_of_Vector2d, (vector2string(obj, 20),)
copy_reg.constructor(construct_vector_of_Vector2d)
copy_reg.pickle(vector_of_Vector2d, reduce_vector_of_Vector2d, construct_vector_of_Vector2d)

def construct_vector_of_Tensor2d(encoded_string):
    return string2vector_of_Tensor2d(encoded_string)
def reduce_vector_of_Tensor2d(obj):
    return construct_vector_of_Tensor2d, (vector2string(obj, 20),)
copy_reg.constructor(construct_vector_of_Tensor2d)
copy_reg.pickle(vector_of_Tensor2d, reduce_vector_of_Tensor2d, construct_vector_of_Tensor2d)

def construct_vector_of_SymTensor2d(encoded_string):
    return string2vector_of_SymTensor2d(encoded_string)
def reduce_vector_of_SymTensor2d(obj):
    return construct_vector_of_SymTensor2d, (vector2string(obj, 20),)
copy_reg.constructor(construct_vector_of_SymTensor2d)
copy_reg.pickle(vector_of_SymTensor2d, reduce_vector_of_SymTensor2d, construct_vector_of_SymTensor2d)

def construct_vector_of_ThirdRankTensor2d(encoded_string):
    return string2vector_of_ThirdRankTensor2d(encoded_string)
def reduce_vector_of_ThirdRankTensor2d(obj):
    return construct_vector_of_ThirdRankTensor2d, (vector2string(obj, 20),)
copy_reg.constructor(construct_vector_of_ThirdRankTensor2d)
copy_reg.pickle(vector_of_ThirdRankTensor2d, reduce_vector_of_ThirdRankTensor2d, construct_vector_of_ThirdRankTensor2d)

def construct_vector_of_Vector3d(encoded_string):
    return string2vector_of_Vector3d(encoded_string)
def reduce_vector_of_Vector3d(obj):
    return construct_vector_of_Vector3d, (vector2string(obj, 20),)
copy_reg.constructor(construct_vector_of_Vector3d)
copy_reg.pickle(vector_of_Vector3d, reduce_vector_of_Vector3d, construct_vector_of_Vector3d)

def construct_vector_of_Tensor3d(encoded_string):
    return string2vector_of_Tensor3d(encoded_string)
def reduce_vector_of_Tensor3d(obj):
    return construct_vector_of_Tensor3d, (vector2string(obj, 20),)
copy_reg.constructor(construct_vector_of_Tensor3d)
copy_reg.pickle(vector_of_Tensor3d, reduce_vector_of_Tensor3d, construct_vector_of_Tensor3d)

def construct_vector_of_SymTensor3d(encoded_string):
    return string2vector_of_SymTensor3d(encoded_string)
def reduce_vector_of_SymTensor3d(obj):
    return construct_vector_of_SymTensor3d, (vector2string(obj, 20),)
copy_reg.constructor(construct_vector_of_SymTensor3d)
copy_reg.pickle(vector_of_SymTensor3d, reduce_vector_of_SymTensor3d, construct_vector_of_SymTensor3d)

def construct_vector_of_ThirdRankTensor3d(encoded_string):
    return string2vector_of_ThirdRankTensor3d(encoded_string)
def reduce_vector_of_ThirdRankTensor3d(obj):
    return construct_vector_of_ThirdRankTensor3d, (vector2string(obj, 20),)
copy_reg.constructor(construct_vector_of_ThirdRankTensor3d)
copy_reg.pickle(vector_of_ThirdRankTensor3d, reduce_vector_of_ThirdRankTensor3d, construct_vector_of_ThirdRankTensor3d)

#------------------------------------------------------------------------------
# std::pairs
#------------------------------------------------------------------------------
# Same as above
## pair_template = """
## def construct_pair_%(value_type1)s_%(value_type2)s(encoded_string1,
##                                                    encoded_string2):
##     return pair_%(value_type1)s_%(value_type2)s(pickle.loads(encoded_string1),
##                                                 pickle.loads(encoded_string2))
## def reduce_pair_%(value_type1)s_%(value_type2)s(obj):
##     return construct_pair_%(value_type1)s_%(value_type2)s, (pickle.dumps(obj.first), pickle.dumps(obj.second))
## copy_reg.pickle(pair_%(value_type1)s_%(value_type2)s, reduce_pair_%(value_type1)s_%(value_type2)s, construct_pair_%(value_type1)s_%(value_type2)s)
## """

## pair_value_types = [("double", "double"),
##                     ("double", "string"),
##                     ("Vector1d", "Vector1d"),
##                     ("Vector2d", "Vector2d"),
##                     ("Vector3d", "Vector3d")]

## for value_type1, value_type2 in pair_value_types:
##     exec(pair_template % {"value_type1": value_type1,
##                           "value_type2": value_type2})

def construct_pair_double_double(encoded_string1,
                                                   encoded_string2):
    return pair_double_double(pickle.loads(encoded_string1),
                                                pickle.loads(encoded_string2))
def reduce_pair_double_double(obj):
    return construct_pair_double_double, (pickle.dumps(obj.first), pickle.dumps(obj.second))
copy_reg.pickle(pair_double_double, reduce_pair_double_double, construct_pair_double_double)

def construct_pair_double_string(encoded_string1,
                                                   encoded_string2):
    return pair_double_string(pickle.loads(encoded_string1),
                                                pickle.loads(encoded_string2))
def reduce_pair_double_string(obj):
    return construct_pair_double_string, (pickle.dumps(obj.first), pickle.dumps(obj.second))
copy_reg.pickle(pair_double_string, reduce_pair_double_string, construct_pair_double_string)

def construct_pair_Vector1d_Vector1d(encoded_string1,
                                                   encoded_string2):
    return pair_Vector1d_Vector1d(pickle.loads(encoded_string1),
                                                pickle.loads(encoded_string2))
def reduce_pair_Vector1d_Vector1d(obj):
    return construct_pair_Vector1d_Vector1d, (pickle.dumps(obj.first), pickle.dumps(obj.second))
copy_reg.pickle(pair_Vector1d_Vector1d, reduce_pair_Vector1d_Vector1d, construct_pair_Vector1d_Vector1d)

def construct_pair_Vector2d_Vector2d(encoded_string1,
                                                   encoded_string2):
    return pair_Vector2d_Vector2d(pickle.loads(encoded_string1),
                                                pickle.loads(encoded_string2))
def reduce_pair_Vector2d_Vector2d(obj):
    return construct_pair_Vector2d_Vector2d, (pickle.dumps(obj.first), pickle.dumps(obj.second))
copy_reg.pickle(pair_Vector2d_Vector2d, reduce_pair_Vector2d_Vector2d, construct_pair_Vector2d_Vector2d)

def construct_pair_Vector3d_Vector3d(encoded_string1,
                                                   encoded_string2):
    return pair_Vector3d_Vector3d(pickle.loads(encoded_string1),
                                                pickle.loads(encoded_string2))
def reduce_pair_Vector3d_Vector3d(obj):
    return construct_pair_Vector3d_Vector3d, (pickle.dumps(obj.first), pickle.dumps(obj.second))
copy_reg.pickle(pair_Vector3d_Vector3d, reduce_pair_Vector3d_Vector3d, construct_pair_Vector3d_Vector3d)
