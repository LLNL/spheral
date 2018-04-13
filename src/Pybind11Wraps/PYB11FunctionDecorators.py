#-------------------------------------------------------------------------------
# Decorators for functions.
#-------------------------------------------------------------------------------
from functools import wraps as PYB11wraps    # Add PYB11 to screen out in generation
import decorator                             # To preserve wrapped functions args

#-------------------------------------------------------------------------------
# cppname
#-------------------------------------------------------------------------------
@decorator.decorator
def PYB11cppname(f, arg):
    @PYB11wraps(f)
    def wrapper(*args, **kwargs):
        return f(*args, **kwargs)
    wrapper.PYB11cppname = arg
    return wrapper

#-------------------------------------------------------------------------------
# pyname
#-------------------------------------------------------------------------------
@decorator.decorator
def PYB11cppname(f):
    @PYB11wraps(f)
    def wrapper(*args, **kwargs):
        return f(*args, **kwargs)
    wrapper.PYB11pyname = arg
    return wrapper

