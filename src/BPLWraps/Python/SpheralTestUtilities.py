# SpheralTestUtilities -- small helper functions used in Spheral unit tests.

import sys
from math import *

import loadmpi
mpi, procID, numProcs = loadmpi.loadmpi()

################################################################################
# Get the frame of the caller of a routine
# Based on Pat Millers getCurrentFrame()
################################################################################
def previousFrame():
    try:
        raise None
    except:
        type, value, traceback = sys.exc_info()
    currentFrame = traceback.tb_frame
    previousFrame = currentFrame.f_back
    return previousFrame

################################################################################
# Get the outermost (presumably global) frame.
# Based on Pat Millers getOuterFrame()
################################################################################
def globalFrame():
    currentFrame = previousFrame()
    globalFrame = currentFrame
    while currentFrame:
        globalFrame = currentFrame
        currentFrame = currentFrame.f_back
    return globalFrame

################################################################################
# Provide a standard title text
################################################################################
def title(titleText, lineLength=80):
    fillerText = "-"*((lineLength - len(titleText))/2)
    print fillerText, titleText, fillerText

################################################################################
# Provide a standard title text
################################################################################
def output(cmd, dict=None):
    if not dict:
        frame = globalFrame()
        dict = frame.f_globals
    print cmd, " : ", eval(cmd, dict)

################################################################################
# Functions to help testing neighbor selection.
################################################################################
def findNeighborNodes(r0, radius, nodes):
    r = nodes.positions()
    H = nodes.Hfield()
    result = [i for i in range(nodes.numInternalNodes)
              if (H[i]*(r[i] - r0)).magnitude() <= radius]
    return result

def checkNeighbors(neighborList, neighborList0):
    if len(neighborList) < len(neighborList0):
        return False
    neighborList.sort()
    neighborList0.sort()
    i = 0
    j = 0
    while i < len(neighborList0) and j < len(neighborList):
        if neighborList[j] == neighborList0[i]:
            i += 1
            j += 1
        elif neighborList[j] > neighborList0[i]:
            return False
        else:
            j += 1
    return i == len(neighborList0) 

################################################################################
# Print statistic about the H tensors for a set of NodeLists.
################################################################################
def hstats(nodeSet):
    for nodes in nodeSet:
        hmin, hmax, havg = 1e50, -1e50, 0.0
        hratmin, hratmax, hratavg = 1e50, -1e50, 0.0
        for H in nodes.Hfield().internalValues():
            Heigen = H.eigenValues()
            hmin = min(hmin, 1.0/(Heigen.maxElement()))
            hmax = max(hmax, 1.0/(Heigen.minElement()))
            havg += sum([1.0/x for x in Heigen.elements()])/Heigen.nDimensions
            hrati = Heigen.minElement()/Heigen.maxElement()
            hratmin = min(hratmin, hrati)
            hratmax = max(hratmax, hrati)
            hratavg += hrati
        n = max(1, mpi.allreduce(nodes.numInternalNodes, mpi.SUM))
        hmin = mpi.allreduce(hmin, mpi.MIN)
        hmax = mpi.allreduce(hmax, mpi.MAX)
        havg = mpi.allreduce(havg, mpi.SUM)/n
        hratmin = mpi.allreduce(hratmin, mpi.MIN)
        hratmax = mpi.allreduce(hratmax, mpi.MAX)
        hratavg = mpi.allreduce(hratavg, mpi.SUM)/n
        print nodes.name(), hmin, hmax, havg, hratmin, hratmax, hratavg

################################################################################
# Grab the node set values (internal, ghost, or all) for FieldLists.
################################################################################
def internalValues(fieldList):
    result = []
    for f in fieldList.fields():
        result.extend(f.internalValues())
    return result

def ghostValues(fieldList):
    result = []
    for f in fieldList.fields():
        result.extend(f.ghostValues())
    return result

def allValues(fieldList):
    result = []
    for f in fieldList.fields():
        result.extend(f.allValues())
    return result

################################################################################
# Fuzzy comparisons.
################################################################################
def fuzzyEqual(lhs, rhs,
               fuzz = 1.0e-5):
    return abs(lhs - rhs)/max(1.0, abs(lhs) + abs(rhs)) < fuzz;

def fuzzyLessThanOrEqual(lhs, rhs,
                         fuzz = 1.0e-5):
    return lhs < rhs or fuzzyEqual(lhs, rhs, fuzz);

def fuzzyGreaterThanOrEqual(lhs, rhs,
                            fuzz = 1.0e-5):
    return lhs > rhs or fuzzyEqual(lhs, rhs, fuzz)

def distinctlyLessThan(lhs, rhs,
                       fuzz = 1.0e-5):
    return lhs < rhs and not fuzzyEqual(lhs, rhs, fuzz)

def distinctlyGreaterThan(lhs, rhs,
                          fuzz = 1.0e-5):
    return lhs > rhs and not fuzzyEqual(lhs, rhs, fuzz)
