# A Spheral utility to help comparing two column oriented files of numbers to some precision.
# Lines starting with "#" are skipped.

import numpy as np

def filearraycmp(filename1, filename2,
                 rtol = 1.0e-5,
                 atol = 1.0e-8):
    array1 = np.loadtxt(filename1)
    array2 = np.loadtxt(filename2)
    #diff = np.absolute(array1 - array2)/np.clip(0.5*(np.absolute(array1) + np.absolute(array2)), atol, 1e100)
    return (np.absolute(array1 - array2) <= atol + rtol*0.5*(np.absolute(array1) + np.absolute(array2))).min()
