#-------------------------------------------------------------------------------
# Find the desired root for a function to the requested accuracy using simple
# bisection.
#-------------------------------------------------------------------------------
def bisectFunction(functor, x1, x2, ytarget, xaccuracy):
    y1 = functor(x1)
    y2 = functor(x2)
    assert (y2 - y1)*(y2 - ytarget) >= 0.0
    while abs(x2 - x1) > xaccuracy:
        xmid = 0.5*(x1 + x2)
        ymid = functor(xmid)
        dymid = ymid - ytarget
        dytar = ytarget - y1
        if dymid*dytar > 0.0:
            x2 = xmid
        else:
            x1 = xmid
            y1 = ymid
    return 0.5*(x1 + x2)
